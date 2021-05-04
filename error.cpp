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
#include "simon_functions.h"


void complexSectors(CoordGrid &gr,  std::vector< int > &activeId, std::vector< int > *sectorTC)
{

  std::vector< GridNode > &Ingrid = gr.m_grid;  
  for(unsigned int n = 0; n < activeId.size(); ++n) {
    int id = activeId[n];
    int idx = gr.Find(id);
    GridNode &curNode = Ingrid[idx];
    int sector = curNode.m_Sector;
    int n_neighbors =0;
    for  ( int i = 0; i < curNode.m_neighbors.size(); i++){
      int neigh_ID = curNode.m_neighbors[i];
      if(neigh_ID < START_VIRTUAL_ID)
	n_neighbors++;	  
    }

    if (n_neighbors > 4 && std::find(sectorTC->begin(), sectorTC->end(),sector)==sectorTC->end() ){
      sectorTC->push_back(sector);
    }
        
  }
}

void complexTracks(CoordGrid &gr,  std::vector< MCTrackObject* >  const *MCTracks, std::vector< int > *idComplex)
{
  std::vector< GridNode > &Ingrid = gr.m_grid;  

  for(size_t j = 0; j < MCTracks->size(); ++j) {
    // Current MC track
    MCTrackObject const *MCtrack = MCTracks->at(j);
    // Stt component of the current MC track.
    std::vector<int> MCSttComp(MCtrack->m_STT_Component);
    for(size_t k = 0; k < MCSttComp.size(); k++){
      int curId = MCSttComp[k];
      int curIdx = gr.Find(curId);
      GridNode &curNode =  Ingrid[curIdx];
      // dbgtrkerror("%d", curId);
      for (size_t i = j+1; i <  MCTracks->size(); i++){
	MCTrackObject const *MCtrack2 = MCTracks->at(i);
	std::vector<int> MCSttComp2(MCtrack2->m_STT_Component);
	for(size_t l =0;l< MCSttComp2.size(); l++){
	  int testId = MCSttComp2[l];
	  if(testId == curId){
	    //    dbgtrkerror("One ID in track %d belongs to the different track %d", MCtrack->m_trackID,MCtrack2->m_trackID);
	    if(std::find(idComplex->begin(), idComplex->end(),MCtrack->m_trackID)==idComplex->end()){
	      //  info("Adding track %d to list of complex", MCtrack->m_trackID);
	      idComplex->push_back(MCtrack->m_trackID);
	    }
	    
	    if(std::find(idComplex->begin(), idComplex->end(),MCtrack2->m_trackID)==idComplex->end()){
	      //  info("Adding track %d to list of complex", MCtrack2->m_trackID);
	      idComplex->push_back(MCtrack2->m_trackID);
	    }
	    break;
	  }
	  if(curNode.IsNeighboring(testId)){
	    //   dbgtrkerror("One ID in track %d is neighbor to the different track %d", MCtrack->m_trackID,MCtrack2->m_trackID);
	    if(std::find(idComplex->begin(), idComplex->end(),MCtrack->m_trackID)==idComplex->end()){
	      // info("Adding track %d to list of complex", MCtrack->m_trackID);
	      idComplex->push_back(MCtrack->m_trackID);
	    }
	    
	    if(std::find(idComplex->begin(), idComplex->end(),MCtrack2->m_trackID)==idComplex->end()){
	      //  info("Adding track %d to list of complex", MCtrack2->m_trackID);
	      idComplex->push_back(MCtrack2->m_trackID);
	    }
	    break;
	  }
	}
      }
      
      if(k > 0 && !curNode.IsNeighboring(MCSttComp[k-1])){
	//	dbgtrkerror(" ID %d in track %d is not neighbor to %d in the same track, complex?", curId, MCtrack->m_trackID, MCSttComp[k-1]);
	if(std::find(idComplex->begin(), idComplex->end(),MCtrack->m_trackID)==idComplex->end()){
	  //	  info("Adding track %d to list of complex", MCtrack->m_trackID);
	  idComplex->push_back(MCtrack->m_trackID);
	}
	break;
      }
    }
  }
  return;
}



MCMatchingError* MatchMCTracksWithConnectedComponents(std::vector< MCTrackObject* >  const *MCTracks,
						      std::vector< std::set<int>* >  const *connectedComp)
{
  if( (MCTracks == 0) || (connectedComp->size() == 0) ) {
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
  
  int nMC = 0;
  // Determine the total area of the all tracks in the current event.
  for(size_t i = 0; i < MCTracks->size(); ++i) {
    //        error("1");

    MCTrackObject const *MCtrack = MCTracks->at(i);
    std::set<int> MCSttComp((MCtrack->m_STT_Component).begin(), (MCtrack->m_STT_Component).end());

    //  error("%ld", MCSttComp.size());

    if(MCSttComp.size()>5){
      TotalArea += (MCtrack->m_STT_Component).size(); //ERROROROROROOR CORRECT AFT CHAP
      nMC++;
    }
  }

  //  error("nMC %d", nMC);
  // FORALL elements in "connectedComp"(T_j) find R_k in "MCTracks"
  // such that (T_j intersection R_k) is maximum.
  std::vector<int>::iterator it;
  std::set<int>::iterator compIt;
  // Component loop(T_j)
  for(size_t i = 0; i < connectedComp->size(); ++i) {
    float matchValue = 0;//std::numeric_limits<float>::min();
    float unionValue = 0;
    float mc_length = 0;

    // Current connected component to analyse
    std::set<int> const *Cur_comp = connectedComp->at(i);
    std::vector<int> Cur_Comp_list;
    for(compIt = Cur_comp->begin(); compIt != Cur_comp->end(); ++compIt) {
      int id = *compIt;
      // Only Real stt tubes, MC does not know anything of the
      // virtuals
      if( (id < START_VIRTUAL_ID ) && (id >= 0) ) {
	Cur_Comp_list.push_back(id);
      }
    }
    if(Cur_Comp_list.size() <= 5)
      continue;
    std::sort(Cur_Comp_list.begin(), Cur_Comp_list.end());
    // Find which MC track has the largest overlap (R_k)
    int idMatchedMC;
    for(size_t j = 0; j < MCTracks->size(); ++j) {
      // Current MC track
      MCTrackObject const *MCtrack = MCTracks->at(j);
      // Stt component of the current MC track.
      std::vector<int> MCSttComp(MCtrack->m_STT_Component);
      std::sort(MCSttComp.begin(), MCSttComp.end());
      // Intersection result output
      std::vector<int> IntersectionList( (Cur_Comp_list.size() + MCSttComp.size()), 0 );

      // Determine the sequence intersection. Only real tubes.
      it = std::set_intersection( Cur_Comp_list.begin(), Cur_Comp_list.end(),
				  MCSttComp.begin(), MCSttComp.end(),
				  IntersectionList.begin());

      // Resize. it points to the position of the last common element
      IntersectionList.resize(it - IntersectionList.begin());
      if(static_cast<float>(IntersectionList.size()) > matchValue) {
	matchValue = IntersectionList.size();
	// A_k
	mc_length = static_cast<float>(MCSttComp.size());

	std::vector<int> UnionList( (Cur_Comp_list.size() + MCSttComp.size()), 0 );

	it=std::set_union ( Cur_Comp_list.begin(), Cur_Comp_list.end(),
				  MCSttComp.begin(), MCSttComp.end(),
			     UnionList.begin());
	UnionList.resize(it - UnionList.begin());
	unionValue = static_cast<float>(UnionList.size());
	idMatchedMC  =  j;
      }
    }// MC tracks loop

    //  dbgtrkerror("Match with MC %d, MatchValue = %f, mc_length = %f",idMatchedMC, matchValue, mc_length);
    
    // We have found ( (T_j intersection R_k) is maximum )
    Error_underMerge +=  ( ( ( mc_length - matchValue ) * matchValue )/ mc_length );
    // Normalized Under segment
    Error_underMergeNorm += ( ( ( mc_length - matchValue ) * matchValue )/ mc_length ) / TotalArea;


    // (alpha_j = Cur_Comp_list.size())
    Error_overMerge  += (static_cast<float>(Cur_Comp_list.size()) - matchValue);

    // Normalized over segment
    Error_overMergeNorm += (static_cast<float>(Cur_Comp_list.size()) - matchValue) / TotalArea;

    //   dbgtrkerror("Undermerge track = %f, over_mergetrack = %f", Error_underMergeNorm,  Error_overMergeNorm);
   
    float Jacard = matchValue / unionValue;

  }// Components loop
  // compute total error for the current set and normalize by the
  // total area of the actual instance (MC Truth)
  TotalError = (Error_underMerge * Error_underMerge) +
               (Error_overMerge  * Error_overMerge);
  // Normalize by total Area of MC.
  TotalError = sqrt(TotalError)/TotalArea;
  
  // Sub terms are already normalised
  TotalErrorNorm = (Error_underMergeNorm * Error_underMergeNorm) +
                   (Error_overMergeNorm  * Error_overMergeNorm);
  TotalErrorNorm = sqrt(TotalErrorNorm);
  
  dbgtrkerror("Error_underMerge = %f, Error_overMerge = %f, TotalError = %f", Error_underMergeNorm, Error_overMergeNorm, TotalErrorNorm);

  if(TotalErrorNorm > 0.7)
    error("Error_underMerge = %f, Error_overMerge = %f, TotalError = %f", Error_underMergeNorm, Error_overMergeNorm, TotalErrorNorm);
  // Prepare output
  MCMatchingError *outError = new MCMatchingError();
  outError->Error_underMerge = Error_underMerge;
  outError->Error_overMerge  = Error_overMerge;
  outError->TotalError       = TotalError;
  // Normalised version
  outError->Error_underMergeNorm = Error_underMergeNorm;
  outError->Error_overMergeNorm  = Error_overMergeNorm;
  outError->TotalErrorNorm       = TotalErrorNorm;
  outError->NumberOfMCTracks     = nMC;
  outError->NumberOfTracklets    = connectedComp->size();

  // Return
  return outError;
}
//__________________ END MatchMCTracksWithConnectedComponents _______________



MCMatchingError* MatchComplexMCTracks(CoordGrid &gr,
				      std::vector< MCTrackObject* >  const *MCTracks,
				      std::vector< std::set<int>* >  const *connectedComp,
				      std::vector<int >  &idComplex )
{
  if( (MCTracks == 0) || (connectedComp == 0) ) {
    error("[COMPLEX] One of the input parameters for matching is empty.");
    return 0;
  }
  std::vector< GridNode > &Ingrid = gr.m_grid;

  // Compute the total MC area (true area) for the current event.
  float TotalArea = 0;
  float Error_underMerge = 0;
  float Error_overMerge  = 0;
  float TotalError = 0;
  // Per error type normalized values
  float Error_underMergeNorm = 0;
  float Error_overMergeNorm  = 0;
  float TotalErrorNorm = 0;
  
  float mc_length = 0;
  bool *include = (bool *) calloc(MCTracks->size(), sizeof(bool));
  // std::vector< bool > include;

  // Determine the total area of the all tracks in the current event.
  for(size_t i = 0; i < MCTracks->size(); ++i) {
    MCTrackObject const *MCtrack = MCTracks->at(i);
    // Stt component of the current MC track.
    std::vector<int> MCSttComp(MCtrack->m_STT_Component);

    if(std::find(idComplex.begin(), idComplex.end(), MCtrack->m_trackID) != idComplex.end()){
      include[i] = true;
      TotalArea += (MCtrack->m_STT_Component).size();
      dbgtrkerror("We will compute an error for the MC track %d with size %d", MCtrack->m_trackID, (MCtrack->m_STT_Component).size());
    }
  }
  // FORALL elements in "connectedComp"(T_j) find R_k in "MCTracks"
  // such that (T_j intersection R_k) is maximum.
  std::vector<int>::iterator it;
  std::set<int>::iterator compIt;
  // Component loop(T_j)
  for(size_t i = 0; i < connectedComp->size(); ++i) {
    float matchValue = std::numeric_limits<float>::min();
    // Current connected component to analyse
    std::set<int> const *Cur_comp = connectedComp->at(i);
    std::vector<int> Cur_Comp_list;
    //   bool include = false;
    for(compIt = Cur_comp->begin(); compIt != Cur_comp->end(); ++compIt) {
      int id = *compIt;
      // Only Real stt tubes, MC does not know anything of the
      // virtuals
      if( (id < START_VIRTUAL_ID ) && (id >= 0) ) 
	Cur_Comp_list.push_back(id);
	/*	int idx = gr.Find(id);
	GridNode* node =&Ingrid[idx];
	int sector = node->m_Sector;
	if(!include && std::find(sectorList.begin(), sectorList.end(), sector) != sectorList.end())
	  include = true;
	  }*/
    }
  //    if (!include) continue;*/
    std::sort(Cur_Comp_list.begin(), Cur_Comp_list.end());
    // Find which MC track has the largest overlap (R_k)
    int id = -1;
    for(size_t j = 0; j < MCTracks->size(); ++j) {
      // Current MC track
      if(!include[j]) continue;
      
      MCTrackObject const *MCtrack = MCTracks->at(j);
      // Stt component of the current MC track.
      std::vector<int> MCSttComp(MCtrack->m_STT_Component);
      std::sort(MCSttComp.begin(), MCSttComp.end());
      // Intersection result output
      std::vector<int> IntersectionList( (Cur_Comp_list.size() + MCSttComp.size()), 0 );

      // Determine the sequence intersection. Only real tubes.
      it = std::set_intersection( Cur_Comp_list.begin(), Cur_Comp_list.end(),
				  MCSttComp.begin(), MCSttComp.end(),
				  IntersectionList.begin());

      // Resize. it points to the position of the last common element
      IntersectionList.resize(it - IntersectionList.begin());
      if(static_cast<float>(IntersectionList.size()) > matchValue) {
	matchValue = IntersectionList.size();
	// A_k
	mc_length = static_cast<float>(MCSttComp.size());
	id = MCtrack->m_trackID;
      }
    }// MC tracks loop

    if(id == -1) continue;
    dbgtrkerror("CC %d, best MatchValue = %f with MC track %d, mc_length = %f", i,  matchValue, id, mc_length);

    
    // We have found ( (T_j intersection R_k) is maximum )
    Error_underMerge +=  ( ( ( mc_length - matchValue ) * matchValue )/ mc_length );
    // Normalized Under segment
    Error_underMergeNorm += ( ( ( mc_length - matchValue ) * matchValue )/ mc_length ) / TotalArea;

    // (alpha_j = Cur_Comp_list.size())
    Error_overMerge  += (static_cast<float>(Cur_Comp_list.size()) - matchValue);

    // Normalized over segment
    Error_overMergeNorm += (static_cast<float>(Cur_Comp_list.size()) - matchValue) / TotalArea;

  }// Components loop
  // compute total error for the current set and normalize by the
  // total area of the actual instance (MC Truth)
  TotalError = (Error_underMerge * Error_underMerge) +
               (Error_overMerge  * Error_overMerge);
  // Normalize by total Area of MC.
  TotalError = sqrt(TotalError)/TotalArea;
  
  // Sub terms are already normalised
  TotalErrorNorm = (Error_underMergeNorm * Error_underMergeNorm) +
                   (Error_overMergeNorm  * Error_overMergeNorm);
  TotalErrorNorm = sqrt(TotalErrorNorm);
  
   dbgtrkerror("[COMPLEX] Error_underMerge = %f, Error_overMerge = %f, TotalError = %f", Error_underMergeNorm, Error_overMergeNorm, TotalErrorNorm);

  // Prepare output
  MCMatchingError *outError = new MCMatchingError();
  outError->Error_underMerge = Error_underMerge;
  outError->Error_overMerge  = Error_overMerge;
  outError->TotalError       = TotalError;
  // Normalised version
  outError->Error_underMergeNorm = Error_underMergeNorm;
  outError->Error_overMergeNorm  = Error_overMergeNorm;
  outError->TotalErrorNorm       = TotalErrorNorm;
  // Return
  free(include);
  return outError;
}
//___

float  calculateAngle(
    float x1, float y1,
    float x2, float y2, 
    float x3, float y3)
{
    // Find direction ratio of line AB
    float ABx = x1 - x2;
    float ABy = y1 - y2;
  
    // Find direction ratio of line BC
    float BCx = x3 - x2;
    float BCy = y3 - y2;
  
    // Find the dotProduct
    // of lines AB & BC
    double dotProduct
        = ABx * BCx
      + ABy * BCy;
  
    // Find magnitude of
    // line AB and BC
    double magnitudeAB
        = ABx * ABx
      + ABy * ABy;
    double magnitudeBC
        = BCx * BCx
      + BCy * BCy;
  
    // Find the cosine of
    // the angle formed
    // by line AB and BC
    double angle = dotProduct;
    angle /= sqrt(
        magnitudeAB * magnitudeBC);
  
    // Find angle in radian
    angle = (angle * 180) / 3.14;
  
    // Print the angle
    //  cout << abs(angle) << endl;
    return (float) angle;
}




//___________________________ MatchPerTrackWithMCTracks ____________________
/* Function to evaluate error per track. MC is used as ground truth. */
std::vector< MCMatchingError* >* MatchPerTrackWithMCTracks(CoordGrid const &hitMap,
                                                           std::vector < MCTrackObject* > const *MCTracks,
                                                           std::vector < PathCandidate* > const *tracklets,
							   std::vector<int >  &idComplex,
							   std::vector<int >  &matchedId)
{
  std::vector< GridNode > Ingrid = hitMap.m_grid;  

  if( (MCTracks == 0) || (tracklets->size() == 0) ) {
    error("One of the input parameters for matching is empty.");
    return 0;
  }
  // Create output parameter
  std::vector< MCMatchingError* >* outPutPar = 0;
  outPutPar = new std::vector< MCMatchingError* >();
  if( outPutPar == 0) {
    error("Could not allocate memory for output list in MatchPerTrackWithMCTracks.");
    delete outPutPar;
    return 0;
  }
  
  //dbgtrkerror("Input MC has %d members and components has %d",MCTracks->size(), tracklets->size());
  // dbgtrkerror("first MC has %d elements and last has %d", ((MCTracks->at(0))->m_STT_Component).size(),(MCTracks->at(MCTracks->size()-1)->m_STT_Component).size() );

  std::vector<int> matchedTrackletIndices;
  std::vector<int> matchedNodes;

  std::vector<int>::iterator FindIndexIt;
  // For all MC-tracklets find the best match in connected componets
  // list.(R_k interset T_j) maximum
  std::vector<int>::iterator it;
  std::set<int>::iterator compIt;
  // MC tracks loop
  for(size_t i = 0; i < MCTracks->size(); ++i) {
    int complex = 0;
    int found50 = 0;
    // Current MC track(R_k)
    MCTrackObject const *MCtrack = MCTracks->at(i);
    // Stt component of the current MC track.
    std::set<int> MCSttComp((MCtrack->m_STT_Component).begin(), (MCtrack->m_STT_Component).end());
    std::vector<point3D> MC3DPt((MCtrack->m_pointSTTCoordList).begin(), (MCtrack->m_pointSTTCoordList).end());
    std::vector<float> isochrone((MCtrack->m_STT_Isochrone).begin(), (MCtrack->m_STT_Isochrone).end());
    std::vector<point3D> MCmomentum((MCtrack->m_STT_Momentum).begin(), (MCtrack->m_STT_Momentum).end());

    // if(MCSttComp.size() < 6)
    //   continue;
    // Not empty MC-track (How is this possible??, empty MC-Tracks.)
    std::vector<int> MCSttCompVect(MCtrack->m_STT_Component);

    // if(MCSttComp.size() < 6)
    //   continue;
    // Not empty MC-track (How is this possible??, empty MC-Tracks.)
      // dbgtrkerror("Size track %d is %lu", i, MCSttComp.size());
      //for(size_t p = 0; p < MCSttCompVect.size(); p++)
	
      //	dbgtrkerror("Contain id  %d", MCSttCompVect[p]);


    int zdetbool = 0;

    //error("%ld", MCSttComp.size());

    if( (MCSttComp.size() > 5) ) {
      if(std::find(idComplex.begin(), idComplex.end(), MCtrack->m_trackID) != idComplex.end())
	complex = 1;
      std::vector<point3D> allmatch, anchorsPts, CCpoints, MCpoints;
      srand(time(NULL)); // Seed the time

      for(compIt = MCSttComp.begin(); compIt != MCSttComp.end(); ++compIt) {
	int id = *compIt;
	int idx = hitMap.Find(id);
	GridNode &node = Ingrid[idx];
	if(node.m_type == GridNode::STT_TYPE_PARA)
	  zdetbool |= 1;
	else if(node.m_type == GridNode::STT_TYPE_SKEW)
	  zdetbool |= 2;

	//	float randomX = ((float) rand()) / (float) RAND_MAX/2.;
	//float randomY = ((float) rand()) / (float) RAND_MAX/2.;

	point3D pt(node.m_x,node.m_y,node.m_z);
       	MCpoints.push_back(pt);
	
      }
      /* for(int ii = 0; ii < MC3DPt.size(); ii++){
	float randomX = ((float) rand()) / (float) RAND_MAX - 0.5;
	float randomY = ((float) rand()) / (float) RAND_MAX - 0.5;
	//	point3D pt(MC3DPt[ii].m_x+randomX,MC3DPt[ii].m_y+randomY,MC3DPt[ii].m_z);
	//	MCpoints.push_back(pt);
	}*/

      //  dbgtrkerror("Contain id  %d", MCSttComp[p]);
      int    matchTrackIndex = -1;
      float  match_length = 0;
      float  mc_length = (float) MCSttComp.size();
      
      float  matchValue = std::numeric_limits<float>::min();
      float  unionValue = std::numeric_limits<float>::min();

      if(std::find(matchedId.begin(), matchedId.end(), MCtrack->m_trackID) != matchedId.end()){
	auto it = find(matchedId.begin(),matchedId.end(),  MCtrack->m_trackID);
	int index = it - matchedId.begin();
	matchTrackIndex = index;
      }
      //  dbgtrkerror("Index found %d", matchTrackIndex);
      
      if(matchTrackIndex >= 0) {
        /* Met of zonder teruglegging */
        // Mark as assigned(to avoid multiple assignements)
        ////matchedTrackletIndices.push_back(matchTrackIndex);
        
        std::set<int> const *bestMatchComponent = (tracklets->at(matchTrackIndex))->m_memberIdSet;
	std::vector<int> const *bestMemberList = (tracklets->at(matchTrackIndex))->m_memberList;

	std::vector<GridNode>  *anchors = &(tracklets->at(matchTrackIndex))->m_anchors;
	//printf("Size anchors %ld \n", anchors->size());
	std::vector<double> const &x = tracklets->at(matchTrackIndex)->m_x;//connectedComp->at(j);
	std::vector<double> const &y = tracklets->at(matchTrackIndex)->m_y;//connectedComp->at(j);
	std::vector<double> const &z = tracklets->at(matchTrackIndex)->m_z;//connectedComp->at(j);

        std::vector<int> bestMatchlist;
        std::vector<float> r_iso;

	int ii =0;
	//	dbgtrkerror("Firts id %d, last id %d",*bestMatchComponent->begin(), *bestMatchComponent->end() );
	/*for(set<int>::iterator p = bestMatchComponent->begin(); p != bestMatchComponent->end(); p++){
	  int element = *p;

	  dbgtrkerror("Best track contains id  %d", element);
	  }*/

	

	
	

        //for(compIt = bestMatchComponent->begin(); compIt != bestMatchComponent->end(); ++compIt) {
	for(size_t pp = 0; pp< bestMemberList->size(); pp++){
          int id = bestMemberList->at(pp);//*compIt;
          // Only Real stt tubes, MC does not know anything of the
          // virtuals
	  point3D pt(x[pp],y[pp],z[pp]);
	  allmatch.push_back(pt);
	  int idx = hitMap.Find(id);
	  GridNode &node = Ingrid[idx];
	  
          if( (id < START_VIRTUAL_ID ) && (id >= 0) ) {
            bestMatchlist.push_back(id);
	    auto it = std::find(MCSttCompVect.begin(), MCSttCompVect.end(), node.m_detID);
	    if (it != MCSttCompVect.end() && node.m_type!=GridNode::STT_TYPE_SKEW)
	      {
		int index = std::distance(MCSttCompVect.begin(), it);
		r_iso.push_back(isochrone[index]);
	      } else
	      r_iso.push_back(0);
	    //  point3D pt2(node.m_x,node.m_y,node.m_z);
	    //  CCpoints.push_back(pt2);
	    // dbgtrkerror(" id  %d should be in both", id);

	  } else
	    r_iso.push_back(0);
	  //	  ii++;
        }
        std::sort(bestMatchlist.begin(), bestMatchlist.end());

	for(size_t ll  = 0; ll <  anchors->size(); ll ++){
	  GridNode &d = anchors->at(ll);
	  //  float randomX = 0.;//((float) rand()) / (float) RAND_MAX/2.;
	  //  float randomY = 0.;//((float) rand()) / (float) RAND_MAX/2.;
	  float randomX = 0;//((float) rand()) / (float) RAND_MAX - 0.5;
	  float randomY = 0;//((float) rand()) / (float) RAND_MAX - 0.5;
	  point3D pt(d.m_xDet+randomX,d.m_yDet+randomY,d.m_z_Det);
	  //  printf("%f, %f, %f\n", d.m_xDet,d.m_yDet,d.m_z_Det);
	  anchorsPts.push_back(pt);
	}

	std::vector<int> IntersectionList( (bestMatchlist.size() + MCSttComp.size()), 0 );
	// Determine the sequence intersection. Only real tubes.
	it = std::set_intersection(bestMatchlist.begin(),bestMatchlist.end(),
				   MCSttComp.begin(), MCSttComp.end(),
				   IntersectionList.begin());
	// Resize. "it" points to the position of the last common element
	IntersectionList.resize(it - IntersectionList.begin());

	matchValue = static_cast<float>(IntersectionList.size());
	// Area of matched tracklet (a_j)
	match_length = static_cast<float>(bestMatchlist.size());

	//	dbgtrkerror(" Match 
	// Index of tracklet with largest intersection
	// matchTrackIndex = j;
	//	matchedNodes = IntersectionList;
	if (matchValue > MCSttComp.size()/2)
	  found50 = 1;
	    
	std::vector<int> UnionList( (bestMatchlist.size() + MCSttComp.size()), 0 );
	// Determine the sequence intersection. Only real tubes.
	it = std::set_union( bestMatchlist.begin(), bestMatchlist.end(),
			     MCSttComp.begin(), MCSttComp.end(),
			     UnionList.begin());
	// Resize. "it" points to the position of the last common element
	UnionList.resize(it - UnionList.begin());
	unionValue = static_cast<float>(UnionList.size());


	//	dbgtrkerror("Match Inter %f, CC length %f, MC length %f", matchValue, match_length, mc_length);
          
	
        // Determine differences
        // All x in MC and not in Tracklet (R_k - T_j)
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

	float TP = matchValue;
	float FP =  static_cast<float>(CurCompdiffMC.size());
	float FN =  static_cast<float>(MCdiffCurComp.size());


	std::vector<double> bestX;
	std::vector<double> bestY;
	std::vector<double> bestZ;
	for(size_t p = 0; p < IntersectionList.size(); p++){
	  //  if(matchValue/unionValue > 0.3){
	  int id = IntersectionList[p];
	  auto it = find((tracklets->at(matchTrackIndex))->m_memberList->begin(),(tracklets->at(matchTrackIndex))->m_memberList->end(),  id);
	  int index = it - (tracklets->at(matchTrackIndex))->m_memberList->begin();
	  bestX.push_back(x[index]);
	  bestY.push_back(y[index]);
	  bestZ.push_back(z[index]);
	}

	//	dbgtrkerror("Size pts is %lu", bestX.size());

	
	double disx =0, disy =0, disz = 0,  mindisx,mindisy, mindisz;
	std::vector<float> alldisx;
	std::vector<float> alldisy;
	std::vector<float> alldisz;
	float alldisxx = 0, alldisyy = 0,alldiszz = 0;
	for(size_t k = 0; k < bestX.size(); k++){
	  double mindis = std::numeric_limits<double>::max();
	  //  dbgtrkerror("%f, %f, %f", bestX[k],bestY[k],bestZ[k]);
	  for(size_t l =0; l < MC3DPt.size(); l++){
	    point3D Mcpt = MC3DPt[l];
	    double curdis = sqrt(pow(bestX[k]-Mcpt.m_x,2)+pow(bestY[k]-Mcpt.m_y,2));
	    if(curdis < mindis) {
	      // dbgtrkerror("Best is %f, %f, %f", Mcpt.m_x,Mcpt.m_y,Mcpt.m_z);
	      mindisx =fabs(bestX[k]-Mcpt.m_x);
	      mindisy =fabs(bestY[k]-Mcpt.m_y);
	      mindisz =fabs(bestZ[k]-Mcpt.m_z);
	      alldisxx = (float) (bestX[k]-Mcpt.m_x);///Mcpt.m_x;
	      alldisyy = (float) (bestY[k]-Mcpt.m_y);///Mcpt.m_y;
	      alldiszz = (float) (bestZ[k]-Mcpt.m_z);///Mcpt.m_z;
	      /*if(fabs(alldisyy) >4){
		error("Mc pt %f", Mcpt.m_y);
		error("Best pt %f, %f, %f", bestY[k]);
		}*/
	      mindis = curdis;
	    }
	  }
	  disx += mindisx;
	  disy += mindisy;
	  disz += zdetbool == 3 ? mindisz: 200;
	  /* if(fabs(alldisyy) >4){
	    error("WHATTTTT? ?? ? ?? ? ");
	    //error("Mc pt %f", Mcpt.m_y);
	    // error("Best pt %f, %f, %f", bestY[k]);
	    }*/
	  alldisx.push_back(alldisxx);
	  alldisy.push_back(alldisyy);
	  if(zdetbool == 3)
	    alldisz.push_back(alldiszz);
	  else
	    alldisz.push_back(200);

	}
	disx/= (double)bestX.size();
	disy/= (double)bestX.size();
	disz/= (double)bestX.size();

	float jacardsingle = matchValue/unionValue;
	float f1score = 2*TP/(2*TP+FP+FN);
	//dbgtrkerror("Xdis %lf, Ydis %lf, Zdis %lf", disx, disy, disz);
	//	dbgtrkerror("%d, %d", MCdiffCurComp.size(), CurCompdiffMC.size());
        // Determine curvature parameters for tracklet and MC
        /*______ Curvature for MC points _____*/
        /* HIER was je mee bezig. Het ziet er OK uit. Andere functie
           om MC data te fitten is gemaakt. */
        //MCtrack->print();// FIXME hier Heb je veranderd(FIXME, 24 nov)
        //std::vector<int> MCSttCompList(MCSttComp.begin(), MCSttComp.end());
        CurvatureParameters mcCurvPars, mcCurvPars2, mcCurvSTT;
        /* De functie hieronder niet voor MC gebruiken. Hij gebruikt
         * de in de grid beschikbare coordinaten en niet de MC punten.*/
        //ComputeCurvatureForListOfnodes(hitMap, MCSttCompList, mcCurvPars);

        //int CircleFit( std::vector<point3D> const &points, CurvatureParameters &curvature);
	std::vector<point3D>  AllSttpoints = MCtrack->m_pointSTTCoordList;
	//	printf("NUMBER PTS IN LIST %d\n", AllSttpoints.size());
	/*	if(MCtrack->m_pointSTTCoordList.size() > mc_length){
	  AllSttpoints.clear();
	  for(int kk = 0; kk < mc_length; kk++){
	    AllSttpoints.insert( AllSttpoints.begin(),MCtrack->m_pointSTTCoordList[kk]);
	  }
	  }*/
	//printf("NUMBER PTS IN LIST %d\n", AllSttpoints.size());
	 mcCurvPars.m_a = -1;
	 mcCurvPars.m_b = -1;
        //std::vector<point3D> const &AllSttpoints = MCtrack->m_pointSTTCoordList;
	 //  int blabla1 = CircleFit(AllSttpoints, mcCurvPars);
	 fit_circle(AllSttpoints, mcCurvPars);
	 fit_circle(MCpoints, mcCurvSTT);

		//	dbgtrkerror("Number of it MC %d", blabla1);

        /*________ Determine curvature for bet match tracklet _____*/
        CurvatureParameters componentCurvPars,componentCurvPars2 ;
        // FIXME: Maybe better to include virtuals. Better gradual ....
	//	ComputeCurvatureForListOfnodes(hitMap, anchorsPts, componentCurvPars);
	componentCurvPars.m_a = mcCurvPars.m_a;
	componentCurvPars.m_b = mcCurvPars.m_b;
	/*int blabla2 = CircleFit(allmatch, componentCurvPars); //allmatch

	if (blabla2 == 4999999){
	  componentCurvPars.m_a = -1;//mcCurvPars.m_a;
	  componentCurvPars.m_b = -1;//mcCurvPars.m_b;
	  blabla2 = CircleFit(allmatch, componentCurvPars);
	  }*/
	//	dbgtrkerror("Number of it CC %d", blabla2);
	//	if(blabla <= 1000)
	//	  error("PROBLEM");
       	fit_circle(anchorsPts, componentCurvPars);
	//float allp = 0;
	std::vector<float> allp;
	for(int in = 0; in < 100; in++){
	  CCpoints.clear();
	  // ii = 0;
	  // for(compIt = bestMatchComponent->begin(); compIt != bestMatchComponent->end(); ++compIt) {
	  for(int kk = 0; kk < r_iso.size();kk++) {
	    // int id = *compIt;
	    // Only Real stt tubes, MC does not know anything of the
	    // virtuals
	    //int idx = hitMap.Find(id);
	    // GridNode &node = Ingrid[idx];
	    // std::vector<int> MCSttCompVect(MCtrack->m_STT_Component);
	    float rx =0, ry =0;
	    //   //rx = ((float) rand()) / (float) RAND_MAX - 0.5;
	    // ry = ((float) rand()) / (float) RAND_MAX - 0.5;
	    // if(r_iso[kk] != 0){// (id < START_VIRTUAL_ID ) && (id >= 0) ) {
	      //auto it = std::find(MCSttCompVect.begin(), MCSttCompVect.end(), node.m_detID);
	      // if (it != MCSttCompVect.end()) && node.m_type!=GridNode::STT_TYPE_SKEW
	      //	{
	      //	  int index = std::distance(MCSttCompVect.begin(), it);
	      //  float r_iso = isochrone[index];

	    if(r_iso[kk] >0){
	      //float distToNext = kk != r_iso.size()-1? sqrt(pow(x[kk+1]-x[kk],2)+pow(y[kk+1]-y[kk],2)):
	      //	sqrt(pow(x[kk]-x[kk-1],2)+pow(y[kk]-y[kk-1],2));
	       float ri = r_iso[kk];
	      // float distMax = sqrt(pow(distToNext,2)+pow(ri,2));
	      // error("%f, %f,  DistNext %f, DistMax %f",x[kk],y[kk],distToNext, distMax);
	      rx = ((float) rand()) / (float) RAND_MAX * ri * 2 -ri;
	      ry =sqrt(pow(ri,2) - pow(rx,2));
	      float randNeg = ((float) rand()) / (float) RAND_MAX - 0.5;
	      // ry = ((float) rand()) / (float) RAND_MAX * ry * 2 -ry;
	      ry = randNeg < 0? -1*ry:ry;
	      //	printf("%f,%f, %f\n", r_iso, rx, ry);
	      //	}
	      /* if( kk != 0 && kk != r_iso.size()-1){

		double maxAng = returnAngle(x[kk-1], x[kk],x[kk+1],y[kk-1], y[kk], y[kk+1]);
		//	float curDist =kk != r_iso.size()-1? sqrt(pow(x[kk+1]-(x[kk]+rx),2)+pow(y[kk+1]-(y[kk]+ry),2)):/
		//	  sqrt(pow(-x[kk-1]+(x[kk]+rx),2)+pow(-y[kk-1]+(y[kk]+ry),2));
		double curAng = returnAngle(x[kk-1], x[kk], x[kk]+rx, y[kk-1], y[kk], y[kk]+ry);

		bool cod = curAng <= MAX(maxAng, 0) && curAng >= MIN(maxAng, 0) || fabs(maxAng) < 170 || fabs(maxAng) > 10? true:false;
		int test = 0;
		while(!cod && test<100){
		  rx = ((float) rand()) / (float) RAND_MAX * ri * 2 -ri;
		  ry =sqrt(pow(ri,2) - pow(rx,2));
		  float randNeg = ((float) rand()) / (float) RAND_MAX - 0.5;
		  // ry = ((float) rand()) / (float) RAND_MAX * ry * 2 -ry;
		  ry = randNeg < 0? -1*ry:ry;
		  //  curDist =kk != r_iso.size()-1? sqrt(pow(x[kk+1]-(x[kk]+rx),2)+pow(y[kk+1]-(y[kk]+ry),2)):
		  //   sqrt(pow(-x[kk-1]+(x[kk]+rx),2)+pow(-y[kk-1]+(y[kk]+ry),2));
		  curAng = returnAngle(x[kk-1], x[kk], x[kk]+rx, y[kk-1], y[kk], y[kk]+ry);
		  cod = curAng <= MAX(maxAng, 0) && curAng >= MIN(maxAng, 0)? true:false;
		  test++;
		}
		
		//	error("Ang max %lf, cur %lf", maxAng,curAng);

		}*/
	      //   error("New pos %f, %f, ",x[kk]+rx,y[kk]+ry);
	
	    }

	    point3D pt(x[kk]+rx,y[kk]+ry,z[kk]);
	    CCpoints.push_back(pt);
	    // }
	    /*  if( (id < START_VIRTUAL_ID ) && (id >= 0) ) {
		bestMatchlist.push_back(id);
		point3D pt2(node.m_x,node.m_y,node.m_z);
		CCpoints.push_back(pt2);
		}*/
	      //ii++;
	  }
	  fit_circle(CCpoints, componentCurvPars2);
	  // if(0.3*2.*1./componentCurvPars2.m_r*0.01 < 20){
	  allp.push_back(componentCurvPars2.m_r);
	  // error("FLY CC rec pt %f", 0.3*2.*1./componentCurvPars2.m_r*0.01);//0.3*2.*1./componentCurvPars2.m_r*0.01

	    //  } else
	//    in--;
	}
	sort( allp.begin(), allp.end() );
	/*	if(f1score == 1 && mc_length >= 20 && 0.3*2.*1./componentCurvPars.m_r*0.01>1){
	  error(" MCtrack %d, MC pt %f", i,0.3*2.*1./mcCurvPars.m_r*0.01);
	  error(" CC rec pt %f, %f, %f", allp[allp.size()/2],allp[84]-allp[allp.size()/2], allp[allp.size()/2]-allp[16]);//0.3*2.*1./componentCurvPars2.m_r*0.01
	error("rec pt %f, rec STT pt %f \n",0.3*2.*1./componentCurvPars.m_r*0.01);
	}*/

	//	error(" CC rec pt %f, %f, %f", allp[allp.size()/2],allp[84]-allp[allp.size()/2], allp[allp.size()/2]-allp[16]);//0.3*2.*1./componentCurvPars2.m_r*0.01
	if(f1score == 1 && mc_length >= 20 && 0.3*2.*1./componentCurvPars.m_r*0.01>2){
	error(" MCtrack %d, MC pt %f, real %f", i,0.3*2.*1./mcCurvPars.m_r*0.01,sqrt(pow(MCmomentum[0].m_x,2) +pow(MCmomentum[0].m_y,2))*2);
	 error("rec pt %f",0.3*2.*1./componentCurvPars.m_r*0.01);
	 error(" CC rec pt %f, %f, %f\n", 0.3*2.*0.01/allp[allp.size()/2],0.3*2.*0.01/allp[16]-0.3*2.*0.01/allp[allp.size()/2],0.3*2.*0.01/ allp[allp.size()/2]-0.3*2.*0.01/allp[84]);//0.3*2.*1./componentCurvPars2.m_r*0.01
	 error("%d", bestMemberList->at(0));
	}
        // Create error object for current sub-path
        MCMatchingError *erroObject = new MCMatchingError();
	erroObject->matchIndex = matchTrackIndex;
        erroObject->Jacardsingle  = jacardsingle;// Not Computed here

	erroObject->disX = disx;
	erroObject->disY = disy;
	erroObject->disZ = disz;
	erroObject->alldisx = alldisx;
	erroObject->alldisy = alldisy;
	erroObject->alldisz = alldisz;
	erroObject->Complex = complex;
        erroObject->BestMatchMCLength  = mc_length ;//MCtrack->m_pointSTTCoordList.size();
        erroObject->CurrentTrackLength = bestMatchlist.size();
        // Set diffs
        erroObject->MCMinCurrentLength = MCdiffCurComp.size();
        erroObject->CurrentMinMCLength = CurCompdiffMC.size();
        // Per track values according to definitions in the book.
        /* {A_k - (T_j ^ R_k)} * (T_j ^ R_k) / A_k */
        erroObject->Error_underMerge = (static_cast<float>(MCSttComp.size()) - matchValue)/unionValue;// Not Computed here
        erroObject->Error_overMerge  = (static_cast<float>(bestMatchlist.size()) - matchValue)/unionValue;// Not Computed here
	//	dbgtrkerror("MCtrack %d, undermerge %f, overmerge %f", i, erroObject->Error_underMerge, erroObject->Error_overMerge );



        erroObject->Jacardaverage  = f1score;// Not Computed here
	
	erroObject->isNotmatched  =  0;
	//dbgtrkerror("MCtrack %d, jacard %f, f1 = %f", i, jacardsingle, f1score);

        // Store curvature parameters. MC

	erroObject->MC_px = (double) MCmomentum[0].m_x;
	erroObject->MC_py = (double) MCmomentum[0].m_y;
	erroObject->MC_pz = (double) MCmomentum[0].m_z;
        erroObject->MC_a  = mcCurvPars.m_a;
        erroObject->MC_b  = mcCurvPars.m_b;
	//	error("%f", sqrt(pow(MCmomentum[0].m_x,2) +pow(MCmomentum[0].m_y,2)));
        erroObject->MC_r = mcCurvPars.m_r;// is 1/r
        erroObject->MC_E =  mcCurvSTT.m_r;//blabla1;
	//	dbgtrkerror("MCtrack %d, MC 1/r %f", i, mcCurvPars.m_r);

        // Matched tracklet
        erroObject->tr_a = allp[84]-allp[allp.size()/2];//componentCurvPars.m_a;
        erroObject->tr_b = allp[allp.size()/2]-allp[16];//componentCurvPars.m_b;
	//,, )

        erroObject->tr_r = componentCurvPars.m_r;// is 1
	erroObject->tr_theta = (tracklets->at(matchTrackIndex))->m_scattAngle;// is 1

	//	if(f1score >= 0.99 && (0.3*2.*1./mcCurvPars.m_r*0.01 > 2 || 0.3*2.*1./componentCurvPars.m_r*0.01 > 2)){
	// error("MCtrack %d, MC pt %f, rec pt %f", i, 0.3*2.*1./mcCurvPars.m_r*0.01, 0.3*2*1./componentCurvPars.m_r*0.01);
	  //  error("Number of it MC %d", blabla1);
	  //  error("Number of it CC %d", blabla2);
	  // error("B 1.5 M2 MCtrack %d, MC pt %f, rec pt %f", i, 0.3*2.*1./mcCurvPars2.m_r*0.01, 0.3*2.*1./componentCurvPars2.m_r*0.01);
	  //  error("Number of it MC 2 %d", mcCurvPars2.m_E);
	  // error("Number of it CC 2 %d\n\n", componentCurvPars2.m_E);


	//	if(f1score == 1 && mc_length  >=20){ // && 0.3*2.*1./componentCurvPars.m_r*0.01 > 5
	//  dbgtrkerror(" MCtrack %d length %f, MC R %f +- %f, MCpts R %f +- %f, rec R %f+- %f, rec STT R %f+- %f,", i,mc_length,1./mcCurvPars.m_r,mcCurvPars.m_E,1./mcCurvSTT.m_r,mcCurvSTT.m_E,1./componentCurvPars.m_r,componentCurvPars.m_E,1./componentCurvPars2.m_r,componentCurvPars2.m_E);
	//dbgtrkerror(" MCtrack %d, MC pt %f, MCpts pt %f, rec pt %f, rec STT pt %f \n", i,0.3*2.*1./mcCurvPars.m_r*0.01,0.3*2.*1./mcCurvSTT.m_r*0.01,0.3*2.*1./componentCurvPars.m_r*0.01, 0.3*2.*1./componentCurvPars2.m_r*0.01);
	  //	}
	//  if (0.3*2.*1./componentCurvPars.m_r*0.01 > 5 && f1score  ==1 && mc_length >=20)
	//   error("ERRROR");
	//	dbgtrkerror("B 1.5 M2 MCtrack %d, MC pt %f, rec pt %f", i, 0.3*2.*1./mcCurvPars2.m_r*0.01, 0.3*2.*1./componentCurvPars2.m_r*0.01);
	  erroObject->tr_E =  allp[allp.size()/2];//componentCurvPars2.m_r;//blabla2;
        erroObject->TotalError       = -1;// Not Computed here
        // Add to output list
        outPutPar->push_back(erroObject);
      }
      else{// not (matchTrackIndex >= 0)
        /* There was not a matching tracklet. Totally mis matched.*/
	
        MCMatchingError *noMatchError = new MCMatchingError();
        noMatchError->isNotmatched = 1.0;
        noMatchError->BestMatchMCLength  = MCSttComp.size();
        noMatchError->CurrentTrackLength = 0;
        noMatchError->MCMinCurrentLength = MCSttComp.size();
        noMatchError->CurrentMinMCLength = 0;
        // Add to output list
	MCMatchingError *erroObject = new MCMatchingError();

	erroObject->disX = -1 ;
	erroObject->disY = -1 ;
	erroObject->disZ = -1 ;
	erroObject->Complex = complex;
	erroObject->BestMatchMCLength  = MCSttComp.size();
	//  erroObject->CurrentTrackLength = bestMatchlist.size();
        // Set diffs
	// erroObject->MCMinCurrentLength = MCdiffCurComp.size();
	// erroObject->CurrentMinMCLength = CurCompdiffMC.size();
        // Per track values according to definitions in the book.
        /* {A_k - (T_j ^ R_k)} * (T_j ^ R_k) / A_k */
        erroObject->Error_underMerge = 1;// Not Computed here
        erroObject->Error_overMerge  = 0;// Not Computed here
	//	dbgtrkerror("MCtrack %d, undermerge %f, overmerge %f", i, erroObject->Error_underMerge, erroObject->Error_overMerge );


        erroObject->Jacardsingle  = 0;// Not Computed here
        erroObject->Jacardaverage = 0;// Not Computed here
	
	erroObject->isNotmatched  = 1;// Not Computed here
	
        // Store curvature parameters. MC
        erroObject->MC_a = 0;
        erroObject->MC_b = 0;
        erroObject->MC_r = 0;// is 1/r
        erroObject->MC_E = 0;
        // Matched tracklet
        erroObject->tr_a = 0;
        erroObject->tr_b = 0;
        erroObject->tr_r = 0;// is 1/r
        erroObject->tr_E = 0;
        erroObject->TotalError = -1;// Not Computed here
        outPutPar->push_back(erroObject);
	//	dbgtrkerror("MSIMATCH: MCtrack %d, MC pt %f, rec pt %f", i, 0.3*1.5*1./ erroObject->MC_r*0.01, 0.3*1.5*1./erroObject->tr_r*0.01);

	// printIntList(MCSttComp);
      }
      
      // dbgtrkerror("MC track %d has best match index %d",MCtrack->m_trackID,matchTrackIndex);
               
    }//END if(MCSttComp.size() != 0) 
    //Go to next MC-Track
  }//MC tracks loop R_j
  // Debug info before return
  // error("Size of per track error list is %d", outPutPar->size());
  return outPutPar;
}
//_______________________ END MatchPerTrackWithMCTracks ___________________
// DEBUG functie, mag weg later. Voegt niets toe.


std::vector< int > BestCompIdToMCTracks( std::vector < MCTrackObject* > const *MCTracks,
					 std::vector < PathCandidate* > const *tracklets)
{

  std::vector<int> matchedId;
  for(size_t i = 0;i< tracklets->size(); i++)
    matchedId.push_back(-1);

  std::vector<int> matchedCCId;
 
  if( (MCTracks == 0) || (tracklets->size() == 0) ) {
    error("One of the input parameters for matching is empty.");
    return matchedId;
  }

  using trMatch = std::pair<int, int>;

  //std::vector<std::vector<int>> trMatch(MCTracks->size());
  //std::vector<std::vector<int>> trID(MCTracks->size());

  std::vector<std::vector<trMatch>> D(MCTracks->size());

  //std::vector<int> FindIndexIt;

  std::vector<int>::iterator FindIndexIt;
  // For all MC-tracklets find the best match in connected componets
  // list.(R_k interset T_j) maximum
  std::vector<int>::iterator it;
  std::set<int>::iterator compIt;



  // MC tracks loop
  for(size_t i = 0; i < MCTracks->size(); ++i) {
    // Current MC track(R_k)
    MCTrackObject const *MCtrack = MCTracks->at(i);
    // Stt component of the current MC track.
    std::set<int> MCSttComp((MCtrack->m_STT_Component).begin(), (MCtrack->m_STT_Component).end());
    std::vector<int> MCSttCompVect(MCtrack->m_STT_Component);

    if( MCSttComp.size() > 5 ) {
 
      int    matchTrackIndex = -1;
      float  match_length = 0;
      float  mc_length = (float) MCSttComp.size();

      float  matchValue = 0;//std::numeric_limits<float>::min();

      for(size_t j = 0; j < tracklets->size(); ++j) {
          std::set<int> const *Cur_comp = (tracklets->at(j))->m_memberIdSet;//connectedComp->at(j);

          std::set<int> Cur_Comp_list;

          for(compIt = Cur_comp->begin(); compIt != Cur_comp->end(); ++compIt) {
            int id = *compIt;
            // Only Real stt tubes, MC does not know anything of the
            // virtuals
            if( (id < START_VIRTUAL_ID ) && (id >= 0) ) {
              Cur_Comp_list.insert(id);
            }
          }//
          // Set intersection
          std::vector<int> IntersectionList( (Cur_Comp_list.size() + MCSttComp.size()), 0 );
          // Determine the sequence intersection. Only real tubes.
          it = std::set_intersection( Cur_Comp_list.begin(), Cur_Comp_list.end(),
                                      MCSttComp.begin(), MCSttComp.end(),
                                      IntersectionList.begin());
          // Resize. "it" points to the position of the last common element
          IntersectionList.resize(it - IntersectionList.begin());

	  // dbgtrkerror("%f", static_cast<float>(IntersectionList.size()));
	  D[i].push_back(make_pair(j, (int) IntersectionList.size()));
	  // dbgtrkerror("i %d (MC id %d), CC %d, match %lu", i, MCtrack->m_trackID,(tracklets->at(j))->m_id, IntersectionList.size());

          /* FIXME hier veranderd om exacte match te vinden */
           if(static_cast<float>(IntersectionList.size()) > matchValue) {
	     //if(static_cast<float>(IntersectionList.size()) == MCSttComp.size()) {
            matchValue = static_cast<float>(IntersectionList.size());
            // Area of matched tracklet (a_j)
	    // dbgtrkerror("CC: First id is %d, last is %d", *( Cur_Comp_list.begin()), *( Cur_Comp_list.end()-1));
            match_length = static_cast<float>(Cur_Comp_list.size());
	    
            // Index of tracklet with largest intersection
            matchTrackIndex = j;
    //	    dbgtrkerror("Size %lu, %lu",Cur_Comp_list.size(),MCSttComp.size());
	    //	    dbgtrkerror("%f, %f", matchValue, unionValue);
          }
	  
	   //  }// End was not assigned before
      }// Components loop (T_j)
      // We found the best matching tracklet for current MC track
      int matchIdx= -1;
      sort(D[i].begin(), D[i].end(), sortbysec2);
      for(size_t p = 0; p < tracklets->size(); p++){
	if(std::find(matchedCCId.begin(), matchedCCId.end(), D[i][p].first) == matchedCCId.end() && D[i][p].second != 0 ){
	  matchedId[D[i][p].first] =MCtrack->m_trackID;
	  matchIdx = D[i][p].first;
	  matchedCCId.push_back(D[i][p].first);
	  break;
	}
      }
     
      dbgtrkerror("MC track %d has best match index %d",MCtrack->m_trackID, matchIdx);

      
    }//END if(MCSttComp.size() != 0) 
    //Go to next MC-Track
  }//MC tracks loop R_j
  // Debug info before return
  return matchedId;
}
//____________
