
#include <iostream>
#include <set>
#include <sstream>
#include <algorithm>
#include  <cmath>
#include <fstream>
#include <vector>
#include <gsl/gsl_poly.h>
// Root headers
#include "TFile.h"
#include "TNtuple.h"
#include "TStopwatch.h"
#include "TH1.h"

// Local headers
#include "auxiliaryfunctions.h"
#include "CollectSttMvdPoints.h"
#include "SttMVDEventDataReader.h"
#include "pathopen.h"
#include "hitcoordinate.h"
#include "floodingFilter.h"
#include "utilfunctions.h"
#include "trackObject.h"
#include "logc.h"
#include "queue.h"
//#include "performFilter.h"
#include "pathCandidate.h"
#include "path_queue.h"


// DEBUG AND STORE definitions

#define EXCLUDE_STTSKEWED_PLOT   0
#define EXCLUDE_VIRTUALS_IN_PLOT 0

#define WRITE_CONNECTED_COMPONENTS 1
#define EVALUATE_ERROR 0
#define INCLUDE_MVD_INOUTPUT_TRACK 0
#define PRINT_TRACKS_DEBUG_INFO 1
#define PRINT_DEBUG_INFO_COMP_MATCH 0
#define PRINT_ERROR_EXTRA_DEBUG_INFO 1
#define WRITE_CONNECTED_COMPONENTS_JSON 0
#define WRITE_CM_ASCII_FILE 1


typedef enum {
  DOWN = 1,
  SAME = 2,
  UP = 4
} Direction;



int determineSkewed_XYPlane_perso( CoordGrid &hitMap, GridNode const &VNode,
                            std::vector<int> &ListOfSkewedNodesIndex,
				   std::vector<int> &ListOfVirtualNodesIndex, char *visited)
{
  info("Correcting xy-coordinates of skewed nodes for this instance for node %d", VNode.m_detID);
	    
  // Fetch all graph nodes
  std::vector< GridNode > &Ingrid = hitMap.m_grid;
  // Neighbours of the input node
  std::vector<int> const &ListOfNeighbors = VNode.m_neighbors;
  // Local variables
  int LastVirtualNodeIdx = -10;
  std::pair<float, float> Rad_deg;
  std::vector<int>::iterator FindIt;
  PathQueue SkewedNodeIndexQueue;
  // There are exactly 2 neighbours (different slopes). And there are
  // two types of virtual nodes (Ax_Vi_Sk or Sk_Vi_S).
  assert( ListOfNeighbors.size() == 2);
  int FNode_id    = ListOfNeighbors[0];//First ID
  int FNode_index = hitMap.Find(FNode_id);//First index
  GridNode &First_Neigh = Ingrid[FNode_index];
  int SNode_id    = ListOfNeighbors[1];//Second ID
  int SNode_index = hitMap.Find(SNode_id);//Second Index
  GridNode &Second_Neigh = Ingrid[SNode_index];
  // All neighbours must be active
  assert( First_Neigh.m_active && Second_Neigh.m_active );
  size_t LocalCurrentLayer, LocalNextLayer;
  LocalCurrentLayer = LocalNextLayer = 0;
  //Axial->Virtual->Skewed
  bool OuterToInner;
  if( visited[FNode_id] < 4) {
    SkewedNodeIndexQueue.inQueue(FNode_index);
    LocalCurrentLayer = First_Neigh.m_Layer;
    info("First node is %d, visited is %d, current layer %d", FNode_id,visited[FNode_id], LocalCurrentLayer);
    OuterToInner =  First_Neigh.m_Layer -  Second_Neigh.m_Layer > 0 ? true : false;
  }

  else{ 
    SkewedNodeIndexQueue.inQueue(SNode_index);
    LocalCurrentLayer = Second_Neigh.m_Layer;
    info("Second node is %d, visited is %d, current layer %d", SNode_id,visited[SNode_id], LocalCurrentLayer);
    OuterToInner =  Second_Neigh.m_Layer - First_Neigh.m_Layer  > 0 ? true : false;


  }

  info("We go in %d direction ( first %d, second %d",OuterToInner, First_Neigh.m_Layer, Second_Neigh.m_Layer);
  //LocalNextLayer = LocalCurrentLayer + 1;
  if(OuterToInner){
    if( LocalCurrentLayer > 0) {
      LocalNextLayer = LocalCurrentLayer - 1;
    }
    else{
      LocalNextLayer = 0;
    }
  }
  else{// if( !OuterToInner) {// Inner to outer
    LocalNextLayer = LocalCurrentLayer + 1;
  }
  // Process skewed node index queue
  while( !SkewedNodeIndexQueue.isEmpty() ) {
    int sk_idx = SkewedNodeIndexQueue.popFront();
    GridNode &SK_Node = Ingrid[sk_idx];
    FindIt = std::find(ListOfSkewedNodesIndex.begin(), ListOfSkewedNodesIndex.end(), sk_idx);
    if( FindIt == ListOfSkewedNodesIndex.end() ) {
      ListOfSkewedNodesIndex.push_back(sk_idx);
    }
    info("Starting with %d", SK_Node.m_detID);
    // List of neighbours
    std::vector<int> const &Neighbours = SK_Node.m_neighbors;
    for( size_t l = 0; l < Neighbours.size(); ++l) {
      int sknID  = Neighbours[l];

      int sknIdx = hitMap.Find(sknID);
      GridNode &SK_NeighNode = Ingrid[sknIdx];
            info("First neighbor is %d on layer %d", sknID, SK_NeighNode.m_Layer);

      // Active and skewed
      if(  (SK_NeighNode.m_type == GridNode::STT_TYPE_SKEW) ) {
        // Current or next layer
        if( (SK_NeighNode.m_Layer == LocalCurrentLayer) ||
            (SK_NeighNode.m_Layer == LocalNextLayer) ) {
	  info("Is skewed and on next layer");
          FindIt = std::find(ListOfSkewedNodesIndex.begin(), ListOfSkewedNodesIndex.end(), sknIdx);
          // Was not added before
          if( FindIt == ListOfSkewedNodesIndex.end() ) {
            SkewedNodeIndexQueue.inQueue(sknIdx);
            ListOfSkewedNodesIndex.push_back(sknIdx);
          }// End If not added before
        }// Current or next layer
      }// END if active and skewed.
      // If neighbour node is a virtual node. We can start
      // correcting and computing the orientations.
      else if( (SK_NeighNode.m_type == GridNode::VIRTUAL_NODE) ) {
	info("Is virtual and on next layer");
	FindIt = std::find(ListOfVirtualNodesIndex.begin(), ListOfVirtualNodesIndex.end(), sknIdx);
	if( (FindIt == ListOfVirtualNodesIndex.end()) && (sknID != VNode.m_detID) ) {
	  // Outer -> inner layer
	  if( (OuterToInner) && (SK_NeighNode.m_Layer <= SK_Node.m_Layer)) {
	    ListOfVirtualNodesIndex.push_back(sknIdx);
	    LastVirtualNodeIdx = sknIdx;
	  }
	  // Inner -> outer layer
	  if( (!OuterToInner) && (SK_NeighNode.m_Layer >= SK_Node.m_Layer)) {
	    ListOfVirtualNodesIndex.push_back(sknIdx);
	    LastVirtualNodeIdx = sknIdx;
	  }
	}
      }// END if active and virtual
    }//END neighbour loop
  }//END WHILE (SkewedNodeIndexQueue)
  int lastVirtualID = -10;
  if (LastVirtualNodeIdx >= 0) {
    GridNode &LastVirtual_node = Ingrid[LastVirtualNodeIdx];
    lastVirtualID = LastVirtual_node.m_detID;
    /* A list of all skewed nodes between the two virtual nodes is
       created by now and we can proceed with corrections of the
       xy-coordinates. Determine the dx and dy between current node and
       last virtual node.*/
    float x_diff = fabs(VNode.m_x - LastVirtual_node.m_x);
    float y_diff = fabs(VNode.m_y - LastVirtual_node.m_y);
    x_diff /= static_cast<float>(ListOfSkewedNodesIndex.size()+1);
    y_diff /= static_cast<float>(ListOfSkewedNodesIndex.size()+1);
    if( VNode.m_x > LastVirtual_node.m_x ) {
      x_diff *= -1;
    }
    if( VNode.m_y > LastVirtual_node.m_y ) {
      y_diff *= -1;
    }
    ////__________________ DEBUG PRINTS ______________________________________
    std::cout << "<DEBUG>\tListOfSkewedNodesIndex.size() " << ListOfSkewedNodesIndex.size()
              << " ListOfVirtualNodesIndex " << ListOfVirtualNodesIndex.size()
              << " x_diff = " << x_diff << ", y_diff = " << y_diff
              << '\n';
    std::cout << " Starting from " << VNode.m_detID << "("<< VNode.m_Layer << "): ";
    for(size_t v = 0; v < ListOfVirtualNodesIndex.size(); ++v) {
      int indexof = ListOfVirtualNodesIndex[v];
      GridNode const &Vnd = Ingrid[indexof];
      std::cout << Vnd.m_detID << "("<< Vnd.m_Layer << "), ";
    }
    std::cout << " Last V_id = " << LastVirtual_node.m_detID << std::endl;
    ////__________________ DEBUG PRINTS ______________________________________
    float xInc = VNode.m_x + x_diff;
    float yInc = VNode.m_y + y_diff;
    /* Correct xy-coordinates of the skewed nodes */
    for(size_t m = 0; m < ListOfSkewedNodesIndex.size(); ++m) {
      GridNode &skewedToproc = Ingrid[ListOfSkewedNodesIndex[m]];
      skewedToproc.m_xDet = xInc;
      skewedToproc.m_yDet = yInc;
      xInc += x_diff;
      yInc += y_diff;
    }
  }
  /* Skewed nodes are pre-processed. RETURN*/
  return(lastVirtualID);
}


bool sortNeighbors(CoordGrid *gr, GridNode *currentNode, std::vector<int> *prev, std::vector<int> *same, std::vector<int> *next, std::vector<int> *virt, char *visited, int *dir){

  int curDir = *dir;
  std::vector< GridNode > &Ingrid = gr->m_grid;
  int curLayer = currentNode->m_Layer;
  int curId =  currentNode->m_detID;
  bool cond = true;

  std::vector<int>  ListOfSkewedNodesIndex;
  std::vector<int>  ListOfVirtualNodesIndex;
  for(int i = 0; i < currentNode->m_neighbors.size(); i++){
    int neighId = currentNode->m_neighbors[i];
    int neighIdx = gr->Find(neighId);
    GridNode *neighNode = &Ingrid[neighIdx];
    debug("Node %d has one neig %d", curId, neighId);
    if(neighNode->m_type == GridNode::VIRTUAL_NODE){
      virt->push_back(neighId);
      //   determineSkewed_XYPlane_perso( *gr, *neighNode, ListOfSkewedNodesIndex, ListOfVirtualNodesIndex,visited);
      //error("%d curlaye %d, virtual layer %d", neighId, curLayer, neighNode->m_Layer);
      continue;
      /*  debug("%d is a virtual node, add to list and find next neighbor", neighId);
      neighId    = neighNode->m_neighbors[0] == curId ?  neighNode->m_neighbors[1]: neighNode->m_neighbors[0];
      neighIdx   = gr->Find(neighId);
      neighNode  = &Ingrid[neighIdx];*/
    }	 

    if(visited[neighId] == 0){
      if(neighNode->m_Layer > curLayer){
	//	debug("Node %d has one neigh up %d", curId, neighId);
	next->push_back(neighId);
	curDir |= UP;
	visited[neighId] = 2;
      } else if( neighNode->m_Layer < curLayer) {
	//	debug("Node %d has one neigh down %d", curId, neighId);
	prev->push_back(neighId);
	curDir |= DOWN;
	visited[neighId] = 2;
      } else {
	//	debug("Node %d has one neigh on the same %d", curId, neighId);
	same->push_back(neighId);
	//	curDir |= SAME;
	visited[neighId] = 2;
      }
    } else if(visited[neighId] == 4){
      //   debug("Node %d has already been connected, tricky", neighId);
      cond = false;
      if(!(std::find(same->begin(), same->end(), neighId) != same->end()))
	same->push_back(neighId);

    }// else
     //  debug("Node %d has already been added, tricky", neighId);
  }

  	
  if( curDir > 4 || same->size() > 1){
    info("Too many possibilities, let skip it for now");
    cond = false;
  }
  *dir = curDir;

  return cond;
  
}

void resetLists(char *visited, std::vector<int> *prev, std::vector<int> *same, std::vector<int> *next){
  for (int i = 0; i < same->size(); i++)
    if(visited[same->at(i)]<4)
      visited[same->at(i)] =0;
	      
  for (int i = 0; i < next->size(); i++)
    if(visited[next->at(i)] < 4)
      visited[next->at(i)] =0;

  for (int i = 0; i < prev->size(); i++)
    if(visited[prev->at(i)] < 4)
      visited[prev->at(i)] =0;
}

void removeIdFromNeigh(GridNode *neighNode, std::vector<int> *prevNodes, int curId){
  for(int i = 0; i < prevNodes->size(); i++){
    if(prevNodes->at(i) != neighNode->m_detID){ 		
      (neighNode->m_neighbors).erase(std::remove((neighNode->m_neighbors).begin(), (neighNode->m_neighbors).end(),prevNodes->at(i)), (neighNode->m_neighbors).end());
      // debug("%d removed from  neighbor list of %d ",prevNodes->at(i), neighNode->m_detID);

      if(i == 0)
	neighNode->parent = curId;
    }
  }
}

bool areAdjacent(CoordGrid *gr, std::vector<int> *v){
  int adjacent = 0;
  std::vector< GridNode > &Ingrid = gr->m_grid;

  for (int i = 0; i < v->size(); i++){
    int neighId   = v->at(i);
    int neighIdx  = gr->Find(neighId);
    GridNode *neighNode = &Ingrid[neighIdx];
    //prevNodes.push_back(neighId);
    
    for (int j = i+1; j < v->size (); j++){
      //   debug("Are %d and %d connected?", neighId, v->at(j));
      if(neighId == v->at(j)) error("HOUSTON");
      else if(neighNode->IsNeighboring(v->at(j))){
	//	debug("Yes");
	adjacent++;
      }// else
      //	debug("No");
    }
  }
  
  if(adjacent >= v->size() -1)
    return true;
  else
    return false;
}


double *polyFit(std::vector<double>  x, std::vector<double>  y){
  int i,j,k,n,N;
  n = 2;
  double X[2*n+1];        //Array that will store the values of sigma(xi),sigma(xi^2),sigma(xi^3)....sigma(xi^2n)
  for (i=0;i<2*n+1;i++)
    {
      X[i]=0;
      for (j=0;j<x.size();j++)
	X[i]=X[i]+pow(x[j],i);        //consecutive positions of the array will store N,sigma(xi),sigma(xi^2),sigma(xi^3)....sigma(xi^2n)
    }
  double B[n+1][n+2];
  double *a = (double*)calloc(n+1, sizeof(double));//B is the Normal matrix(augmented) that will store the equations, 'a' is for value of the final coefficients
  for (i=0;i<=n;i++)
    for (j=0;j<=n;j++)
      B[i][j]=X[i+j];            //Build the Normal matrix by storing the corresponding coefficients at the right positions except the last column of the matrix
  double Y[n+1];                    //Array to store the values of sigma(yi),sigma(xi*yi),sigma(xi^2*yi)...sigma(xi^n*yi)
  for (i=0;i<n+1;i++)
    {    
      Y[i]=0;
      for (j=0;j<x.size();j++)
        Y[i]=Y[i]+pow(x[j],i)*y[j];        //consecutive positions will store sigma(yi),sigma(xi*yi),sigma(xi^2*yi)...sigma(xi^n*yi)
    }
  for (i=0;i<=n;i++)
    B[i][n+1]=Y[i];                //load the values of Y as the last column of B(Normal Matrix but augmented)
  n=n+1;                //n is made n+1 because the Gaussian Elimination part below was for n equations, but here n is the degree of polynomial and for n degree we get n+1 equations
  /*  cout<<"\nThe Normal(Augmented Matrix) is as follows:\n";    
  for (i=0;i<n;i++)            //print the Normal-augmented matrix
    {
      for (j=0;j<=n;j++)
	cout<<B[i][j]<<setw(16);
      cout<<"\n";
      }    */
  for (i=0;i<n;i++)                    //From now Gaussian Elimination starts(can be ignored) to solve the set of linear equations (Pivotisation)
    for (k=i+1;k<n;k++)
      if (B[i][i]<B[k][i])
	for (j=0;j<=n;j++)
	  {
	    double temp=B[i][j];
	    B[i][j]=B[k][j];
	    B[k][j]=temp;
	  }
     
  for (i=0;i<n-1;i++)            //loop to perform the gauss elimination
    for (k=i+1;k<n;k++)
      {
	double t=B[k][i]/B[i][i];
	for (j=0;j<=n;j++)
	  B[k][j]=B[k][j]-t*B[i][j];    //make the elements below the pivot elements equal to zero or elimnate the variables
      }
  for (i=n-1;i>=0;i--)                //back-substitution
    {                        //x is an array whose values correspond to the values of x,y,z..
      a[i]=B[i][n];                //make the variable to be calculated equal to the rhs of the last equation
      for (j=0;j<n;j++)
	if (j!=i)            //then subtract all the lhs values except the coefficient of the variable whose value                                   is being calculated
	  a[i]=a[i]-B[i][j]*a[j];
      a[i]=a[i]/B[i][i];            //now finally divide the rhs by the coefficient of the variable to be calculated
    }
   
  cout<<"The fitted Polynomial is given by:\ny=";
  for (i=0;i<n;i++)
    cout<<" + ("<<a[i]<<")"<<"x^"<<i;
  cout<<"\n";
  return a;
}

int returnDirection(double prev, double cur){
  double diff = cur - prev;
  int dir;
  if (diff > 1.)
    dir = 1;
  else if (diff < -1.)
    dir = -1;
  else
    dir = 0;


  return dir;
}


int fitnextId(CoordGrid &gr, std::vector<double> &x, std::vector<double> &y, std::vector<int> &next, int method){
  
  std::vector< GridNode > &Ingrid = gr.m_grid;

  if(method == 0){ // linear

  }  else if(method == 1){ // quadratic
    int xDir = returnDirection(  x[x.size()-2],  x[x.size()-1]);
    int yDir = returnDirection(  y[y.size()-2],  y[y.size()-1]);

    debug("We are going in x dir %d and y dir %d", xDir, yDir);

    double minDist = std::numeric_limits<double>::max();
    int goodId = -1;
    std::vector<double> p;
    p.push_back(0.);
    
    for (int i = 0; i < x.size() - 1; i++){
      double newval = p[i] + sqrt(pow(x[i+1]-x[i],2.) + pow(y[i+1]-y[i],2.));
      p.push_back(newval);
      // debug("%lf, %lf", p[i+1]);

    }


    double *x_coef = polyFit(p, x);
    double *y_coef = polyFit(p, y);

    for (int i = 0; i <next.size(); i++){
      int curId = next[i];
      int curIdx = gr.Find(curId);
      GridNode *node = &Ingrid[curIdx];
      double xdet = (double) node->m_xDet;
      double ydet = (double) node->m_yDet;

      double newx_coef[3] = {x_coef[0] - xdet, x_coef[1], x_coef[2]};
      double newy_coef[3] = {y_coef[0] - ydet, y_coef[1], y_coef[2]};

      double vectortanx[3] = {x_coef[1], 2*x_coef[2], 0.};
      double vectortany[3] = {y_coef[1], 2*y_coef[2], 0.};
	      
      double d[4] = {newx_coef[0]*vectortanx[0] + newy_coef[0]*vectortany[0],
		     newx_coef[0]*vectortanx[1] + newx_coef[1]*vectortanx[0] +
		     newy_coef[0]*vectortany[1] + newy_coef[1]*vectortany[0],
		     newx_coef[0]*vectortanx[2] + newx_coef[1]*vectortanx[1] + newx_coef[2]*vectortanx[0] +
		     newy_coef[0]*vectortany[2] + newy_coef[1]*vectortany[1] + newy_coef[2]*vectortany[0],
		     newx_coef[2]*vectortanx[1] + newy_coef[2]*vectortany[1]};

      double x0[3];
      int nroot = gsl_poly_solve_cubic(d[2]/d[3], d[1]/d[3], d[0]/d[3], x0, x0+1, x0+2);

      //debug("Polynomial coeff : %lf + %lf x + %lf x^2 + %lf x^3", d[0], d[1], d[2], d[3]);
      //	debug("Real roots %d : x1 = %lf x2 =  %lf x3 =  %lf \n\n", nroot, x0[0], x0[1], x0[2]);

      double xIntersect, yIntersect, currDist;
      if(nroot ==1){
	xIntersect = gsl_poly_eval(x_coef, 3, x0[0]);
	yIntersect = gsl_poly_eval(y_coef, 3, x0[0]);
	currDist = sqrt(pow(xIntersect - xdet,2) + pow(yIntersect - ydet,2));
      }

		
      for (int j = 0; j < nroot; j++){
	double newx = gsl_poly_eval(x_coef, 3, x0[j]);
	double newy = gsl_poly_eval(y_coef, 3, x0[j]);
	double newdist = sqrt(pow(newx - xdet,2) + pow(newy - ydet,2));
	if(j == 0){
	  xIntersect = newx;
	  yIntersect = newy;
	  currDist = newdist;
	} else if( currDist > newdist) {
	  xIntersect = newx;
	  yIntersect = newy;
	  currDist = newdist;
	}
      }

      debug("Distance between point %d and fit is %lf", curId, currDist);
      // double disttube = distanceBetweenTube(node, lastNode);
      //debug("Tube Distance to point %d is %lf",curId, disttube);


      int newxDir = returnDirection(  x[x.size()-1], xdet);
      int newyDir = returnDirection(  y[y.size()-1], ydet);

      debug("This node is in x dir %d and y dir %d", newxDir, newyDir);

		
      if(minDist > currDist){
	if(labs(newxDir - xDir) > 1 || labs(newyDir - yDir) > 1){
	  debug("Let's avoid going back if we can");
	} else{
		    
	  minDist = currDist;
		  
	  goodId = curId;
	  //goodNode = ;
	}
      }
	    
		
    } // FOR RIGHT NEIGHBORS
    info("The good id is %d, and is at distance %lf", goodId,minDist);
    return goodId;

  }  

  return -1;
}




//______________________ BEGIN MCTrackPoints _____________________________
std::vector< std::vector < MCTrackObject* >* >*
MCTrackPoints( std::vector < std::vector<HitCoordinate*>* > const &evtData)
{
  info("Extracting MC tracks for %d events", evtData.size());
  //	    << " events.\n";
  // Output Parameter
  std::vector< std::vector < MCTrackObject* >* >* outVar =
    new std::vector< std::vector < MCTrackObject* >* >();

  int numTracks = 0;
  std::vector< int > idtracks;

   for(size_t e = 0; e < evtData.size(); ++e) {
    std::vector<HitCoordinate*> const *Current_Event = evtData[e];
    // numTracks = 0;
    // Find out how many MC tracks are available
    for(size_t i = 0; i < Current_Event->size(); ++i) {
      HitCoordinate const *currentHit = Current_Event->at(i);
      if( currentHit->m_trackID != HIT_EXCLUSION && !(std::find(idtracks.begin(), idtracks.end(), currentHit->m_trackID) != idtracks.end())) {
	printf("%d \n", currentHit->m_trackID);
	//	numTracks = currentHit->m_trackID;
	idtracks.push_back(currentHit->m_trackID);
      }
    }// END current event loop
    numTracks = idtracks.size();
    debug("Event %d contains %d tracks", e,(numTracks));
    
   // Now, we know the number of available tracks for the current
    // event; We can allocate memory.
    std::vector < MCTrackObject* >* evtTracks = new std::vector < MCTrackObject* >();
    for(int j = 0; j < numTracks; ++j) {
      MCTrackObject *trk = new MCTrackObject(); 
      evtTracks->push_back(trk);
    }
    // Fill the list

    int curTrack = 0;
    for(size_t k = 0; k < Current_Event->size(); ++k) {
      HitCoordinate const *currentHit = Current_Event->at(k);
      if(currentHit->m_trackID != HIT_EXCLUSION) {
	std::vector<int>::iterator it = std::find(idtracks.begin(), idtracks.end(), currentHit->m_trackID);
	int index = std::distance(idtracks.begin(), it);
	//	printf("%d, %d \n", currentHit->m_trackID,index);
	int trackPos = std::distance(idtracks.begin(), it); //currentHit->m_trackID;
	point3D spacePoint;
	spacePoint.m_x = currentHit->mx;
	spacePoint.m_y = currentHit->my;
	spacePoint.m_z = currentHit->mz;
	if(currentHit->type == HitCoordinate::STT_TYPE) {
	  ((evtTracks->at(trackPos))->m_pointSTTCoordList).push_back(spacePoint);
	  ((evtTracks->at(trackPos))->m_STT_Component).push_back(currentHit->m_detID);
	}
	else if(currentHit->type == HitCoordinate::MVD_TYPE) {
	  ((evtTracks->at(trackPos))->m_pointMVDCoordList).push_back(spacePoint);
	  ((evtTracks->at(trackPos))->m_MVD_Component).push_back(currentHit->m_detID);
	}
      }// END if not HIT_EXCLUSION
    }// END
    outVar->push_back(evtTracks);
  }// END events loop
  return outVar;
  /*  int realnumTracks = 0;

  std::vector<int> IDtracks ;


  for(size_t e = 0; e < evtData.size(); ++e) {
    std::vector<HitCoordinate*> const *Current_Event = evtData[e];
    // numTracks = 0;
    // Find out how many MC tracks are available
    for(size_t i = 0; i < Current_Event->size(); ++i) {
      HitCoordinate const *currentHit = Current_Event->at(i);
      if (std::find(IDtracks.begin(), IDtracks.end(), currentHit->m_trackID) == IDtracks.end())
	IDtracks.push_back(currentHit->m_trackID);

      if( currentHit->m_trackID > numTracks &&
	  currentHit->m_trackID != HIT_EXCLUSION) {
	numTracks = currentHit->m_trackID;
      }

    }// END current event loop

#if (PRINT_TRACKS_DEBUG_INFO > 0)
    std::cout << "\t<-I-> Event " << e
	      << " Contains " << (IDtracks.size())
	      << " Tracks.\n";
#endif
    // Now, we know the number of available tracks for the current
    // event; We can allocate memory.
    std::vector < MCTrackObject* >* evtTracks = new std::vector < MCTrackObject* >();
    for(int j = 0; j <= numTracks; ++j) {
      MCTrackObject *trk = new MCTrackObject(); 
      evtTracks->push_back(trk);
    }
    // Fill the list
    for(size_t k = 0; k < Current_Event->size(); ++k) {
      HitCoordinate const *currentHit = Current_Event->at(k);
      if(currentHit->m_trackID != HIT_EXCLUSION) {
	int trackPos = currentHit->m_trackID;
	point3D spacePoint;
	spacePoint.m_x = currentHit->mx;
	spacePoint.m_y = currentHit->my;
	spacePoint.m_z = currentHit->mz;
	if(currentHit->type == HitCoordinate::STT_TYPE) {
	  ((evtTracks->at(trackPos))->m_pointSTTCoordList).push_back(spacePoint);
	  ((evtTracks->at(trackPos))->m_STT_Component).push_back(currentHit->m_detID);
	}
	else if(currentHit->type == HitCoordinate::MVD_TYPE) {
	  ((evtTracks->at(trackPos))->m_pointMVDCoordList).push_back(spacePoint);
	  ((evtTracks->at(trackPos))->m_MVD_Component).push_back(currentHit->m_detID);
	}
      }// END if not HIT_EXCLUSION
    }// END
    outVar->push_back(evtTracks);
  }// END events loop
  return outVar;*/
}
//______________________ END MCTrackPoints _______________________________







void floodingFilter(std::string const &OutFileName,int firstEvt, int lastEvt)
{
  TStopwatch timer;

  // Structure to hold the detector data (grid)
  std::vector < GridNode > detNodes;
  // File and structure to hold the output coordinates.
  TFile Out_Put_File(OutFileName.c_str(),"RECREATE","Outputfile Created by performFilter", 9);

  // Collected coordinates for tracks.
  TNtuple coord ("CoordCollected" , "Collected Coordinates in x y plane", "x:y:z:x_Det:y_Det:z_Det");

  // Ntuple to hold Error values for all events available in the
  // current events set. The value is evaluated per image.
  std::string errorParameter = "Error_underMerge:Error_overMerge:TotalError";
  errorParameter += ":Error_underMergeNorm:Error_overMergeNorm:TotalErrorNorm";
  // Create Ntuple to hold parameters.
  TNtuple ErrorNtuple("ErrorEstimate","Segmentation error values", errorParameter.c_str());

  // Second error type. Per track error value. Based on curvature
  // data.
  std::string PerTrkErrPars = "misMatched:BestMatchMCLength:CurrentTrackLength";
  PerTrkErrPars += ":MCMinCurrentLength:CurrentMinMCLength";
  PerTrkErrPars += ":UnderMergeError:OverMergeError:MC_a:MC_b:MC_r:MC_E:tr_a:tr_b:tr_r:tr_E";
  TNtuple ErrorNtuplePerTrack("PerTrackError","Per track values of error", PerTrkErrPars.c_str());
  
  // NTuple to hold the coordinates of all connected components.
  std::string ConnCompPar = "EvtNum:CompNum:tubeId:x:y:z:x_Det:y_Det:z_Det";
  TNtuple ConnectedCoord ("ConnectedComponents", "Connected component Coordinates",
			  ConnCompPar.c_str());
  // Hold number of components per event
  TNtuple ComponentPerEvt ("ComponentPerEvt", "Component per event","evtNum:numComponents");


  /* Read all data directly from sim, digi and parameter files
     (OLD) */
  std::vector < std::vector<HitCoordinate*>* >* Hit_coords = 
    CollectSttMvdPoints(detNodes, Out_Put_File, firstEvt, lastEvt);

  std::vector< std::vector < MCTrackObject* >* > *MC_Tracks = MCTrackPoints(*Hit_coords);
  
  // Write event info to output  
  WriteEventPlotsToFile( (*Hit_coords), Out_Put_File);
  
  // Create an empty grid object
  CoordGrid gr;
  
  // Init Grid for STT detector nodes (fill the map).
  gr.Initialize(detNodes);

  TNtuple Layers("LayerLimits","Layer Limits.","x:y:det_z:z");
  TNtuple Sections("SectionsLimits","Section Limits.","x:y:det_z:z");
  // Isolate Sector and Layer Limits
  isolateSectorAndLayerLimits(gr, Sections, Layers);
  Sections.Write();
  Layers.Write();

  TNtuple* OrigGrid = GridToNtuple(detNodes, "OrigGridCoord");
  OrigGrid->SetMarkerStyle(8);
  OrigGrid->SetMarkerSize(0.2);
  OrigGrid->SetMarkerColor(kBlack);
  OrigGrid->Write();

  info("Fix neighbouring before extension.");
  fixNeighboring(gr);
  std::vector < GridNode > VNodes;
  ///Compute_Add_VirtualNodes_Neigbor2(gr, VNodes);// Old method
  Compute_Add_VirtualNodes_Neigbor(gr, VNodes);

  // Compute_Virtual_InterSector_Nodes(gr, 6,VNodes);
   TNtuple* virtualTubes = GridToNtuple(VNodes, "VirtualNodes");
  virtualTubes->SetMarkerStyle(8);
  virtualTubes->SetMarkerSize(0.2);
  virtualTubes->SetMarkerColor(kMagenta);
  virtualTubes->Write();
      /*
   * Extend the grid with new virtual points and fix the missing
   * neighboring relations
   */
  // Extend the grid by virtual nodes between the layer.


  info("Extending the grid by %d virtual nodes between the layers.", VNodes.size());
  gr.ExtendedGrid(VNodes);
  Fix_InterSector_Nodes(gr, 6);

  info("Fixing neighbouring after extension.");
  fixNeighboring(gr);
  TNtuple* extendedGrid = GridToNtuple(gr.m_grid, "ExtendedGrid");
  extendedGrid->SetMarkerStyle(7);//8
  extendedGrid->SetMarkerSize(0.3);
  extendedGrid->SetMarkerColor(17);//41
  extendedGrid->Write();

  info("Total number of tubes after extension = %d", gr.GetNumNodes());
  // Delete allocated memory
  delete (OrigGrid);
  delete virtualTubes;
  delete extendedGrid;
  unsigned int totalnumEvt = Hit_coords->size();
  // Start the timer.

  if(true){

  std::vector< std::set<int>* >* connectedComp = 0;
  connectedComp = new std::vector< std::set<int>* >();
  
  timer.Start();

  /* Fill the grid with the current hit points and process. Handles
     each event separately.*/
  // Event loop

  for(size_t k = 0; k < Hit_coords->size(); ++k) {
    // Data for the current event
    info("Processing event: %d", k);
    std::vector<HitCoordinate*> const *dd = 0;
    dd = Hit_coords->at(k);
    if(dd) {
      gr.FillGrid(*dd);
    }

    std::vector< GridNode > &Ingrid = gr.m_grid;  
    std::vector< int > activeId;
    std::vector< int > remainingActiveId;


    /* Pushing all active detectors into queue */
    
    int nactive = 0, nactive_queue = 0;
    for(unsigned int n = 0; n < Ingrid.size(); ++n) {
      if( Ingrid[n].m_active){
	int NodeId = Ingrid[n].m_detID;
	if(Ingrid[n].m_type != GridNode::VIRTUAL_NODE  ) {
	  activeId.push_back(NodeId);
	  nactive_queue++;
	}
	nactive++;
      }
      
    }

    info("Found the active detectors");

    /* Keep in neighbors only active ones */
    
    for(unsigned int n = 0; n < nactive_queue; ++n) {
      int curid    = activeId[n];
      int curindex = gr.Find(curid);
      GridNode &current_Node = Ingrid[curindex];
      int n_neighbors = current_Node.m_neighbors.size();
      //  debug("%d, has %d neighbors", current_Node.m_detID, n_neighbors);
      for  ( int i = 0; i < current_Node.m_neighbors.size(); i++){
	int neigh_ID = current_Node.m_neighbors[i];
	//	debug("%d curr", neigh_ID);
	int neigh_index = gr.Find(neigh_ID);
	if(!Ingrid[neigh_index].m_active){
	  //  debug("Removing");
	  (current_Node.m_neighbors).erase((current_Node.m_neighbors).begin()+i);
	  i--;
	}
      }
      n_neighbors = current_Node.m_neighbors.size();
      //    debug("%d, has %d neighbors", current_Node.m_detID, n_neighbors);

    }

    info("Removed non active neighbors");

    
    int candidateId = 0;
    std::vector < PathCandidate* > tracklets;
    char *visited = (char *)calloc(50000, sizeof(int));

    info("First step, let's find the obvious tracks");

    
    for(unsigned int n = 0; n < nactive_queue; ++n) {
      
      std::vector<int> sameLayer;
      std::vector<int> nextLayer;
      std::vector<int> prevLayer;
      std::vector<int> nextVirt;
      std::vector<int> prevNodes;
      std::vector<int> *v;

      int curId = activeId[n];
      int curIdx = gr.Find(curId);
      GridNode *currentNode = &Ingrid[curIdx];
      
      int dir = 0;

      if(visited[curId] == 0 && ((currentNode->m_LayerLimit == 1 && currentNode->m_neighbors.size() <= 2) ||  currentNode->m_neighbors.size() == 1)){	
	int n_neighbors = currentNode->m_neighbors.size();
	int curLayer = currentNode->m_Layer;

	info("\n Starting node %d, neigh %d", curId, n_neighbors);

	int neighId;
	int neighIdx;
	GridNode *neighNode;
	bool cond = true;

	cond = sortNeighbors(&gr, currentNode, &prevLayer, &sameLayer, &nextLayer, &nextVirt, visited,  &dir);
	
	if(cond == false){
	  info("Too many possibilities, let skip it for now");
	  resetLists(visited, &prevLayer, &sameLayer, &nextLayer);
	  continue;
	}
	
	PathCandidate *cand = new PathCandidate();// Create a new candidate
	cand->m_id = candidateId++;// Set id
	cand->m_isValid = true;
	cand->insertNewNode(currentNode, 1);
	cand->m_firstNodeVisited = curId;
	if(abs(currentNode->m_SectorLimit) > 0)
	  cand->m_isOnSectorLimit = true;
	visited[curId] = 4;
	prevNodes.push_back(curId);

	int n_connected = 0;

	while(cond){


	  /* 1 NEIGHBOR */
	  
	  if(n_neighbors == 1){
	    //  info("dir %d", dir);
	    if (dir == UP){
	      v = &nextLayer;
	      info("Going up from node %d to node %d", curId, v->at(0));
	    } else if (dir == DOWN){
	      v = &prevLayer;
	      info("Going down from node %d to node %d", curId, v->at(0));
	    } else {
	      v = &sameLayer;
	      info("Going same level from node %d to node %d", curId, v->at(0));
	    }
	     
	    
	    if(nextVirt.size() > 0){
	      // debug("Handling the virtuals");
	      
	      for (int i = 0; i <nextVirt.size(); i++){
		neighId   = nextVirt[i];
		neighIdx  = gr.Find(neighId);
		neighNode = &Ingrid[neighIdx];
		cand->insertNewNode(neighNode,1);
		visited[neighId] = 4;
		n_connected++;

		removeIdFromNeigh(neighNode, &prevNodes, curId);
	      }
	      
	      // prevNodes.clear();
	      prevNodes.insert(prevNodes.end(),  nextVirt.begin(),  nextVirt.end());
	      nextVirt.clear();
	    }
	      
	    

	    neighId    = v->at(0);
	    neighIdx   = gr.Find(neighId);
	    neighNode  = &Ingrid[neighIdx];
	    cand->insertNewNode(neighNode,1);
	    visited[neighId] = 4;
	    n_connected++;
	    
	    removeIdFromNeigh(neighNode, &prevNodes, curId);

	    curId       = neighId;
	    curIdx      = neighIdx;
	    currentNode = neighNode;
	    curLayer    = currentNode->m_Layer;
	    
	    nextLayer.clear();  sameLayer.clear();   prevLayer.clear();
	    dir = 0;
	    prevNodes.clear();
	    prevNodes.push_back(curId);
	    
	    cond = sortNeighbors(&gr, currentNode, &prevLayer, &sameLayer, &nextLayer, &nextVirt, visited,  &dir);

	    n_neighbors = sameLayer.size() + prevLayer.size() + nextLayer.size();
	    info("%d nodes were connected,  %d found for the next step, cond %d \n", n_connected, n_neighbors, cond);
	    n_connected = 0;

	    
	  }


	  /* 1 Same layer neighbor */

	  else if (sameLayer.size() > 0){

	    v = dir == UP ? &nextLayer: &prevLayer;

	    if(!areAdjacent(&gr, v)){
	      info("Too many possibilities");
	      cond = false;
	      //resetLists(visited, &prevLayer, &sameLayer, &nextLayer);
	      // break;
	    } else { 
	    
	      int candId = sameLayer[0];
	      int candIdx        = gr.Find(candId);
	      GridNode *candNode = &Ingrid[candIdx];
	      curLayer           = candNode->m_Layer;
	      sameLayer.clear();

	      info("Still on the same layer, investigating node %d", candId);


	      bool toadd = true;
	      int cursize = v->size();

	    
	      for(int i = 0; i < candNode->m_neighbors.size(); i++){
		neighId = candNode->m_neighbors[i];
		if (neighId == curId) continue;
	      
		neighIdx = gr.Find(neighId);
		neighNode = &Ingrid[neighIdx];
	     
		if(neighNode->m_type == GridNode::VIRTUAL_NODE){
		  //	nextVirt.push_back(neighId);
		  /* neighId = neighNode->m_neighbors[0] == candId ? neighNode->m_neighbors[1]: neighNode->m_neighbors[0];
		     neighIdx = gr.Find(neighId);
		     neighNode  = &Ingrid[neighIdx];*/
		  continue;
		}
	      
		int this_neigh = 0;
		for (int j  = 0; j < v->size(); j++){
		  debug("Connection between %d and %d ?", neighId, v->at(j));

		  if(v->at(j) == neighId){
		    //	  debug("They are the same");
		    this_neigh = 1;
		  }
		  else if(neighNode->IsNeighboring(v->at(j))){
		    //	  debug("They are neighbors",  neighId, v->at(j));
		    this_neigh = 1;
		  } else {
		    debug("%d and %d are not (obviously) connected", v->at(j), neighId);
		  }
		}
		if(this_neigh == 0){
		  toadd = false;
		  cond = false;
		  break;
		}
	      } // for in neighbors
		
	      if (toadd == true){
		info("All neighbors look good !");

		cand->insertNewNode(candNode,1);
		visited[candId] = 4;
		n_connected++;

		removeIdFromNeigh(candNode, &prevNodes, curId);
	      
		prevNodes.push_back(candId);
		curId       = candId;
		curIdx      = candIdx;
		currentNode = candNode;
		cond = sortNeighbors(&gr, currentNode, &prevLayer, &sameLayer, &nextLayer, &nextVirt, visited,  &dir);

		n_neighbors = sameLayer.size()+prevLayer.size()+nextLayer.size();
		info("%d nodes were connected,  %d found for the next step \n\n", n_connected, n_neighbors);
		n_connected = 0;
	      }
	      else {

		info("Neighbors not connected... \n\n", n_connected, n_neighbors);
		visited[candId] = 0;
		cand->m_rpotNeigh.push_back(candId);
		info("visited curID %d, canId %d", visited[curId], visited[candId]);

		//   resetLists(visited, &prevLayer, &sameLayer, &nextLayer);

	      }
	    }
	  }



	  /* MORE NEIGHBOOOORS */




	  else if (n_neighbors > 1) {

	    //  info("%d", dir);
	    if (dir == UP){	      
	      info("Next nodes are up and we have %d of them", nextLayer.size());
	      v = &nextLayer;
	    }  else if (dir == DOWN){
	      info("NExt nodes are down and we have %d of them", prevLayer.size());
	      v = &prevLayer;
	    } else
	      info("WHAT IS THE DIRECTION NOW?");
	    
	    
	    if(areAdjacent(&gr, v)){
	      info("Adding %d nodes to the CM", v->size());

	      if(nextVirt.size() > 0){
		//	debug("Handling the virtual");
		for (int i = 0; i < nextVirt.size(); i++){
		  neighId   = nextVirt[i];
		  neighIdx  = gr.Find(neighId);
		  neighNode = &Ingrid[neighIdx];
		  cand->insertNewNode(neighNode,1);
		  visited[neighId] = 4;
		  n_connected++;

		  removeIdFromNeigh(neighNode, &prevNodes, curId);		  
		}
		
		//	prevNodes.clear();
		prevNodes.insert(prevNodes.end(),  nextVirt.begin(),  nextVirt.end());
		nextVirt.clear();

	      }
	      
	      std::vector<int> lookneigh;

	      for (int i = 0; i < v->size(); i++){
		neighId = v->at(i);
		neighIdx = gr.Find(neighId);
		neighNode = &Ingrid[neighIdx];
		cand->insertNewNode(neighNode,1);
		visited[neighId] = 4;
		n_connected++;
		lookneigh.push_back(neighId);

		removeIdFromNeigh(neighNode, &prevNodes, curId);		  
		removeIdFromNeigh(neighNode, v, curId);		  

	      }
	      
	      nextLayer.clear();
	      sameLayer.clear();
	      prevLayer.clear();
	      prevNodes.clear();
	      
	      dir = 0;
	     
	      curId       = lookneigh[0];
	      curIdx      = gr.Find(curId);
	      currentNode = &Ingrid[curIdx];
	      curLayer    = currentNode->m_Layer;
	      info("New current node is %d, looking for neighbors of %d nodes", curId, lookneigh.size());
	      
	      n_neighbors = 0;
	      
	      for(int i = 0; i < lookneigh.size(); i++){
		int id         = lookneigh[i];
		int idx        = gr.Find(id);
		GridNode *node = &Ingrid[idx];
		prevNodes.push_back(id);

		cond = sortNeighbors(&gr, node, &prevLayer, &sameLayer, &nextLayer, &nextVirt, visited,  &dir);
	      }

	      n_neighbors = sameLayer.size()+prevLayer.size()+nextLayer.size();
	      info("%d nodes were connected,  %d found for the next step \n", n_connected, n_neighbors);
	      n_connected = 0;
	     

	    } else {
	      info("Some of these nodes are node adjacent", v->size());
	      cond = false;
	    }
	    
	      
	  } else {

	    if(n_neighbors == 0) {
	      info("No more neighbors in sight, are we finished here ? \n\n");
	       if(cand->m_minLayer == 0 && cand->m_maxLayer > 21){
		 info("track goes through all layers, likily finished");
		 cand->m_finished = 1;
	       } else if(abs(currentNode->m_SectorLimit) > 0 || cand->m_isOnSectorLimit){
		 
		 cand->m_isOnSectorLimit= true;
		 info("Track is on sector limit, might have a connection somewhere else");
	       } else {
		 int firstId = cand->m_firstNodeVisited;
		 int firstIdx = gr.Find(firstId);
		 GridNode *firstNode = &Ingrid[firstIdx];
		 if(firstNode->m_LayerLimit == 1 && currentNode->m_LayerLimit == 1) {
		   info("starting and ending on same layer");
		   cand->m_finished = 1;

		 }
	       }
	    } else {
	      info("This track has other neighbors, to correct later");
	    }
	    cond = false;
	    
	  }

	  debug("Cond %d, neigh %d", cond, n_neighbors);
	  
	  if(cond == false){
	    info("We are going out of the boucle");

	    for (int i = 0; i < sameLayer.size(); i++)
	      cand->m_rpotNeigh.push_back(sameLayer[i]);
	      
	    for (int i = 0; i < nextLayer.size(); i++)
	      cand->m_rpotNeigh.push_back(nextLayer[i]);

	    for (int i = 0; i < prevLayer.size(); i++)
	      cand->m_rpotNeigh.push_back(prevLayer[i]);
	    
	    resetLists(visited, &prevLayer, &sameLayer, &nextLayer);
	    //  continue;
	  }

	}// WHILE COND
	
       	info("Pushing cm %d: \n               length is %d, \n               Min layer %d, \n              Max layer %d, \n               IsOnSectorLimit %d, \n               Finished ? %d. ", cand->m_id, cand->m_length,  cand->m_minLayer, cand->m_maxLayer, cand->m_isOnSectorLimit, cand->m_finished);

	tracklets.push_back(cand);

      }
      // elser
      //	remainingActiveId.push_back(curId);
       	     
    }





    
    if(false){
 
    info(" Find remaining obvious tracklets");

    for(unsigned int n = 0; n < activeId.size(); ++n) {
      int curId = activeId[n];
      int oriId = curId;
      int prevId = curId;
      int curIdx = gr.Find(curId);
      GridNode *currentNode = &Ingrid[curIdx];
      GridNode *firstNode = &Ingrid[curIdx];
      std::vector<int> *nextIds;// List of members in a vector(delete me)

      int n_neighbors = currentNode->m_neighbors.size();

      if(n_neighbors == 2 && visited[curId] == 0 && currentNode->m_LayerLimit == 0) {

	int neighId1 = currentNode->m_neighbors[0];
	int neighIdx1 = gr.Find(neighId1);
	GridNode *neighNode1 = &Ingrid[neighIdx1];

	int neighId2 = currentNode->m_neighbors[1];
	int neighIdx2 = gr.Find(neighId2);
	GridNode *neighNode2 = &Ingrid[neighIdx2];
	
	if(!neighNode1->IsNeighboring(neighId2) && visited[neighId1] == 0 && visited[neighId2] == 0){
	  
	  info("%d has two neighbors not connected %d %d, sounds like a good tracklet", curId, neighId1, neighId2);

	
	  PathCandidate *cand = new PathCandidate();// Create a new candIdate
	  cand->m_id = candidateId++;// Set id
	  cand->m_isValid = true;
	  cand->insertNewNode(currentNode,1);
	  visited[curId] = 4;

	  for (int k = 0; k < 2; k++){
	    std::vector<int> sameLayer;
	    std::vector<int> nextLayer;
	    std::vector<int> prevLayer;
	    std::vector<int> nextVirt;
	    std::vector<int> prevNodes;
	    std::vector<int> *v;

	    int dir = 0;
      
	    int curLayer = firstNode->m_Layer;
	    
	    int neighId =  firstNode->m_neighbors[k];
	    int neighIdx = gr.Find(neighId);
	    GridNode *neighNode = &Ingrid[neighIdx];
	    
	    info("Handling  node %d", neighId);

	    cand->insertNewNode(neighNode,k==0? 0:1);
	    visited[neighId] = 4;
	    
	    (neighNode->m_neighbors).erase(std::remove((neighNode->m_neighbors).begin(), (neighNode->m_neighbors).end(), oriId), (neighNode->m_neighbors).end());

	    //  prevId = curId;
	    curId = neighId;
	    curIdx = neighIdx;
	    currentNode = neighNode;
	    curLayer = currentNode->m_Layer;
	    
	    bool cond = true;
	    
	    cond = sortNeighbors(&gr, currentNode, &prevLayer, &sameLayer, &nextLayer, &nextVirt, visited,  &dir);

	    
	    if(cond == false){
	      //     if(k == 0) cand->m_firstNodeVisited = curId;
	      //else cand->m_lastNodeVisited = curId;

	      	if(k == 0){
		  cand->m_firstNodeVisited = curId;
		  for (int i = 0; i < sameLayer.size(); i++)
		    cand->m_lpotNeigh.push_back(sameLayer[i]);
	      
		  for (int i = 0; i < nextLayer.size(); i++)
		    cand->m_lpotNeigh.push_back(nextLayer[i]);

		  for (int i = 0; i < prevLayer.size(); i++)
		    cand->m_lpotNeigh.push_back(prevLayer[i]);
		}
		else {
		  for (int i = 0; i < sameLayer.size(); i++)
		    cand->m_rpotNeigh.push_back(sameLayer[i]);
	      
		  for (int i = 0; i < nextLayer.size(); i++)
		    cand->m_rpotNeigh.push_back(nextLayer[i]);

		  for (int i = 0; i < prevLayer.size(); i++)
		    cand->m_rpotNeigh.push_back(prevLayer[i]);
		  cand->m_lastNodeVisited = curId;
		}
	      // add neighbors
	      resetLists(visited, &prevLayer, &sameLayer, &nextLayer);

	      continue;
	    }

	  

	    prevNodes.push_back(curId);
	    n_neighbors = prevLayer.size() + sameLayer.size() + nextLayer.size();
	    int n_connected = 0;

	    info("Next number of neighbors is %d", n_neighbors);
	    
	    while(cond){


	      /* 1 NEIGHBOR */
	  
	      if(n_neighbors == 1){
		info("dir %d", dir);
		if (dir == UP){
		  v = &nextLayer;
		  info("Going up from node %d to node %d", curId, v->at(0));
		} else if (dir == DOWN){
		  v = &prevLayer;
		  info("Going down from node %d to node %d", curId, v->at(0));
		} else {
		  v = &sameLayer;
		  info("Going same level from node %d to node %d", curId, v->at(0));
		}
	     
	    
		if(nextVirt.size() > 0){
		  debug("Handling the virtuals");
	      
		  for (int i = 0; i <nextVirt.size(); i++){
		    neighId   = nextVirt[i];
		    neighIdx  = gr.Find(neighId);
		    neighNode = &Ingrid[neighIdx];
		    cand->insertNewNode(neighNode,k==0? 0:1);
		    visited[neighId] = 4;
		    n_connected++;

		    removeIdFromNeigh(neighNode, &prevNodes, curId);
		  }
	      
		  //	  prevNodes.clear();
		  prevNodes.insert(prevNodes.end(),  nextVirt.begin(),  nextVirt.end());
		  nextVirt.clear();
		}
	      
	    

		neighId    = v->at(0);
		neighIdx   = gr.Find(neighId);
		neighNode  = &Ingrid[neighIdx];
		cand->insertNewNode(neighNode,k==0? 0:1);
		visited[neighId] = 4;
		n_connected++;
	    
		removeIdFromNeigh(neighNode, &prevNodes, curId);

		curId       = neighId;
		curIdx      = neighIdx;
		currentNode = neighNode;
		curLayer    = currentNode->m_Layer;
	    
		nextLayer.clear();  sameLayer.clear();   prevLayer.clear();
		dir = 0;
		prevNodes.clear();
		prevNodes.push_back(curId);
	    
		cond = sortNeighbors(&gr, currentNode, &prevLayer, &sameLayer, &nextLayer, &nextVirt, visited,  &dir);

		n_neighbors = sameLayer.size() + prevLayer.size() + nextLayer.size();
		info("%d nodes were connected,  %d found for the next step \n", n_connected, n_neighbors);
		n_connected = 0;

	    
	    
	      }


	      /* 1 Same layer neighbor */

	      else if (sameLayer.size() > 0){

		v = dir == UP ? &nextLayer: &prevLayer;
		if(!areAdjacent(&gr, v)){
		  info("Too many possibilities");
		  //resetLists(visited, &prevLayer, &sameLayer, &nextLayer);
		  //break;
		  cond = false;
		} else{
		  int candId = sameLayer[0];
		  int candIdx        = gr.Find(candId);
		  GridNode *candNode = &Ingrid[candIdx];
		  curLayer           = candNode->m_Layer;
		  sameLayer.clear();

		  info("Still on the same layer, investigating node %d", candId);


		  bool toadd = true;
		  int cursize = v->size();

		  for(int i = 0; i < candNode->m_neighbors.size(); i++){
		    neighId = candNode->m_neighbors[i];
		    if (neighId == curId) continue;
	      
		    neighIdx = gr.Find(neighId);
		    neighNode = &Ingrid[neighIdx];
	     
		    if(neighNode->m_type == GridNode::VIRTUAL_NODE){
		      nextVirt.push_back(neighId);
		      /* neighId = neighNode->m_neighbors[0] == candId ? neighNode->m_neighbors[1]: neighNode->m_neighbors[0];
			 neighIdx = gr.Find(neighId);
			 neighNode  = &Ingrid[neighIdx];*/
		      continue;
		    }
	      


		    int this_neigh = 0;
		    for (int j  = 0; j < v->size(); j++){
		      debug("Connection between %d and %d ?", neighId, v->at(j));

		      if(v->at(j) == neighId){
			debug("They are the same");
			this_neigh = 1;
		      }
		      else if(neighNode->IsNeighboring(v->at(j))){
			debug("They are neighbors",  neighId, v->at(j));
			this_neigh = 1;
		      } else {
			debug("%d and %d are not (obviously) connected", v->at(j), neighId);
		      }
		    }
		    if(this_neigh == 0){
		      toadd = false;
		      cond = false;
		      break;
		    }
		  } // for in neighbors
		
		  if (toadd == true){
		    info("All neighbors look good !");

		    cand->insertNewNode(candNode,k==0? 0:1);
		    visited[candId] = 4;
		    n_connected++;

		    removeIdFromNeigh(candNode, &prevNodes, curId);
	      
		    prevNodes.push_back(candId);
		    curId       = candId;
		    curIdx      = candIdx;
		    currentNode = candNode;
		    cond = sortNeighbors(&gr, currentNode, &prevLayer, &sameLayer, &nextLayer, &nextVirt, visited,  &dir);

		    n_neighbors = sameLayer.size()+prevLayer.size()+nextLayer.size();
		    info("%d nodes were connected,  %d found for the next step \n\n", n_connected, n_neighbors);
		    n_connected = 0;
		  }
		  else {

		    info("To many possibilities, we'll see what's next \n\n", n_connected, n_neighbors);

		    visited[candId] = 0;
		    if (k == 0)
		      cand->m_lpotNeigh.push_back(candId);
		    else
		      cand->m_rpotNeigh.push_back(candId);

		    //		  resetLists(visited, &prevLayer, &sameLayer, &nextLayer);
		    //		  if(k == 0) cand->m_firstNodeVisited = curId;
		    //		  else cand->m_lastNodeVisited = curId;	      
		  

		  }
		}
	      }



	      /* MORE NEIGHBOOOORS */




	      else if (n_neighbors > 1) {

		info("%d", dir);
		if (dir == UP){	      
		  info("Next nodes are up and we have %d of them", nextLayer.size());
		  v = &nextLayer;
		}  else if (dir == DOWN){
		  info("NExt nodes are down and we have %d of them", prevLayer.size());
		  v = &prevLayer;
		} else
		  info("WHAT IS THE DIRECTION NOW?");
	    
	    
		if(areAdjacent(&gr, v)){
		  info("Adding %d nodes to the CM", v->size());

		  if(nextVirt.size() > 0){
		    debug("Handling the virtual");
		    for (int i = 0; i < nextVirt.size(); i++){
		      neighId   = nextVirt[i];
		      neighIdx  = gr.Find(neighId);
		      neighNode = &Ingrid[neighIdx];
		      cand->insertNewNode(neighNode,k==0? 0:1);
		      visited[neighId] = 4;
		      n_connected++;

		      removeIdFromNeigh(neighNode, &prevNodes, curId);		  
		    }
		
		    // prevNodes.clear();
		    prevNodes.insert(prevNodes.end(),  nextVirt.begin(),  nextVirt.end());
		    nextVirt.clear();

		  }
	      
		  std::vector<int> lookneigh;

		  for (int i = 0; i < v->size(); i++){
		    neighId = v->at(i);
		    neighIdx = gr.Find(neighId);
		    neighNode = &Ingrid[neighIdx];
		    cand->insertNewNode(neighNode,k==0? 0:1);
		    visited[neighId] = 4;
		    n_connected++;
		    lookneigh.push_back(neighId);

		    removeIdFromNeigh(neighNode, &prevNodes, curId);		  
		    removeIdFromNeigh(neighNode, v, curId);		  

		  }
	      
		  nextLayer.clear();
		  sameLayer.clear();
		  prevLayer.clear();
		  prevNodes.clear();
	      
		  dir = 0;
	     
		  curId       = lookneigh[0];
		  curIdx      = gr.Find(curId);
		  currentNode = &Ingrid[curIdx];
		  curLayer    = currentNode->m_Layer;
		  info("New current node is %d, looking for neighbors of %d nodes", curId, lookneigh.size());
	      
		  n_neighbors = 0;
	      
		  for(int i = 0; i < lookneigh.size(); i++){
		    int id         = lookneigh[i];
		    int idx        = gr.Find(id);
		    GridNode *node = &Ingrid[idx];
		    prevNodes.push_back(id);

		    cond = sortNeighbors(&gr, node, &prevLayer, &sameLayer, &nextLayer, &nextVirt, visited,  &dir);
		  }

		  n_neighbors = sameLayer.size()+prevLayer.size()+nextLayer.size();
		  info("%d nodes were connected,  %d found for the next step \n", n_connected, n_neighbors);
		  n_connected = 0;
	     

		} else {
		  info("Some of these nodes are node adjacent", v->size());
		  cond = false;
		}
	      
	      } else {

		if(n_neighbors == 0) {
		  info("No more neighbors in sight, are we finished here ? \n\n");
		  if(cand->m_minLayer == 0 && cand->m_maxLayer > 21){
		    info("track goes through all layers, likily finished");
		    cand->m_finished = 2;
		  } else if(abs(currentNode->m_SectorLimit) > 0 || cand->m_isOnSectorLimit){
		 
		    cand->m_isOnSectorLimit= true;
		    info("Track is on sector limit, might have a connection somewhere else");
		    cand->m_finished = 1;

		  } else {
		    int firstId = cand->m_firstNodeVisited;
		    int firstIdx = gr.Find(firstId);
		    GridNode *firstNode = &Ingrid[firstIdx];
		    if(firstNode->m_LayerLimit == 1 && currentNode->m_LayerLimit == 1) {
		      info("starting and ending on same layer");
		      cand->m_finished = 2;

		    }
		  }
		} else {
		  info("This track has other neighbors, to correct later");
		}
		cond = false;
	    
	      }

	      if(cond == false){
		info("We are going out of the boucle");
		if(k == 0){
		  cand->m_firstNodeVisited = curId;
		  for (int i = 0; i < sameLayer.size(); i++)
		    cand->m_lpotNeigh.push_back(sameLayer[i]);
	      
		  for (int i = 0; i < nextLayer.size(); i++)
		    cand->m_lpotNeigh.push_back(nextLayer[i]);

		  for (int i = 0; i < prevLayer.size(); i++)
		    cand->m_lpotNeigh.push_back(prevLayer[i]);
		}
		else {
		  for (int i = 0; i < sameLayer.size(); i++)
		    cand->m_rpotNeigh.push_back(sameLayer[i]);
	      
		  for (int i = 0; i < nextLayer.size(); i++)
		    cand->m_rpotNeigh.push_back(nextLayer[i]);

		  for (int i = 0; i < prevLayer.size(); i++)
		    cand->m_rpotNeigh.push_back(prevLayer[i]);
		  cand->m_lastNodeVisited = curId;
		}
		  
		resetLists(visited, &prevLayer, &sameLayer, &nextLayer);

		//continue;
	      }
	      
	    }
	

	  }
       	     
	

	  tracklets.push_back(cand);
	  info("Pushing cm %d: \n               length is %d, \n               Min layer %d, \n              Max layer %d, \n               IsOnSectorLimit %d, \n               Finished ? %d. ", cand->m_id, cand->m_length,  cand->m_minLayer, cand->m_maxLayer, cand->m_isOnSectorLimit, cand->m_finished);

	}

      }
    }


    }






    

    info("Refinement");

    std::sort(tracklets.begin(), tracklets.end(), compareTwoPathsLength); 
    for(unsigned int l = 0; l < tracklets.size(); l++){
      PathCandidate &curCand = *(tracklets[l]);
      info("Length %d", curCand.m_length);
      //  computePathCurvature(gr,curCand);
      //info("%lf", (curCand.m_CurV_par).m_r);
    }

    char *visitedTracks = (char *) calloc(tracklets.size(), sizeof(char));

    if(true){
      
    for(unsigned int l = 0; l < tracklets.size(); l++){
      PathCandidate &curCand = *(tracklets[l]);
      if (curCand.m_finished > 0){
	info("This CM is either finished or has no further close neighbors");
	//check that it has not pot neighbors
	
	debug("size lpotneigh %d, size rpotneigh %d", curCand.m_lpotNeigh.size(), curCand.m_rpotNeigh.size());
      } else {
	int first =  curCand.m_firstNodeVisited;
	int firstIdx = gr.Find(first);
	GridNode *firstNode = &Ingrid[firstIdx];
	int n_neighfirst = firstNode->m_neighbors.size();

	int last =  curCand.m_lastNodeVisited;
	int lastIdx = gr.Find(last);
	GridNode &lastNode = Ingrid[lastIdx];
	int n_neighlast = lastNode.m_neighbors.size();
	info("This cm has first node %d with %d neighbors and last node %d with %d neighbors", first, n_neighfirst, last, n_neighlast);

	info("left neighbors size %d", curCand.m_lpotNeigh.size());
	if(curCand.m_lpotNeigh.size() > 0)
	  for(int i = 0; i  < curCand.m_lpotNeigh.size(); i++)
	    debug("%d", curCand.m_lpotNeigh[i]);

	info("right neighbors size %d", curCand.m_rpotNeigh.size());
	if(curCand.m_rpotNeigh.size() > 0)
	  for(int i = 0; i  < curCand.m_rpotNeigh.size(); i++)
	    debug("%d", curCand.m_rpotNeigh[i]);

	if(curCand.m_length > 6 && (curCand.m_rpotNeigh.size() > 0 || curCand.m_lpotNeigh.size() > 0 )){
	  info("Finding potential continuation");


	  //storing x and y detections in vectors



	  if(curCand.m_rpotNeigh.size() > 0 && visitedTracks[curCand.m_id] == 0){
	    info("Handling right neighbors, %d of them", curCand.m_rpotNeigh.size());
	    std::vector<int> next;
	    next.insert(next.end(),  (curCand.m_rpotNeigh).begin(),  (curCand.m_rpotNeigh).end());

	    bool cond = true;
	    int potCm = -1;
	    
	    std::vector<int>  *trk = curCand.m_memberList;

	    std::vector<double> x;
	    std::vector<double> y;
	    for (int i = trk->size() - MIN(trk->size(), 10) ; i < trk->size() ; i++){
	      int id = trk->at(i);
	      int idx = gr.Find(id);
	      printf("%d \t", id);
	      GridNode &node = Ingrid[idx];
	      // if(node.m_type == GridNode::VIRTUAL_NODE)
		
	      x.push_back(node.m_xDet);
	      y.push_back(node.m_yDet);
	    }
	    printf("\n");

	    // finding previous direction //
	    
	    while (cond){	      
	      GridNode *goodNode;
	      int goodId =fitnextId(gr, x, y, next, 1);

	      if(goodId != -1){
		int goodIdx = gr.Find(goodId);
		goodNode = &Ingrid[goodIdx];
		debug("CHECK GOOOD IF %d", goodId);
		debug("CHECK GOOOD IF %d", goodNode->m_detID);

		if(visited[goodId] == 4){
		  info(", also this node belongs to this cm %d", goodNode->m_cm[0]);

		  if(std::find((curCand.m_toMerge).begin(), (curCand.m_toMerge).end(), goodNode->m_cm[0]) != (curCand.m_toMerge).end()){
		    debug("We were already planning to merge these tracks, so we stop");
		    cond = false;
		  }
		  else if(potCm == -1)
		    potCm = goodNode->m_cm[0];
		  else if( potCm != goodNode->m_cm[0]){
		    debug("We were visiting different cm, weird");
		    potCm = goodNode->m_cm[0];
		  } else {
		    debug("Two nodes in a row belong to the same cm %d, let's stop and merge them later", potCm);
		    cond = false;
		    (curCand.m_toMerge).push_back(potCm);
		    (tracklets[potCm]->m_toMerge).push_back(curCand.m_id);
		    visitedTracks[curCand.m_id] = 1;
		    visitedTracks[potCm] = 1;
		  }
		} else
		  potCm = -1;
	      } else {
		info("No good candidates have been found, stop");
		cond = false;
	      }
	      
	      if(cond){
		// Find next neighbors
		//	int nonVisited = 0;
		//std::vector<int> potTracks;
		next.clear();
		for(int i = 0; i < goodNode->m_neighbors.size(); i++){
		  int neighId = goodNode->m_neighbors[i];
		  debug("Node %d has one neig %d, visited %d", goodId, neighId, visited[neighId]);

		  int neighIdx = gr.Find(neighId);
		  GridNode *neighNode = &Ingrid[neighIdx];
		  if(curCand.isInCandidate(neighId)) continue;

		  if(neighNode->m_type == GridNode::VIRTUAL_NODE)
		    continue;
		  
		  if(visited[neighId] == 4){
		    debug("Node %d has already been connected, let see", neighId);
		    // potTracks.push_back(neighNode->m_cm[0]);
		  }

		    next.push_back(neighId);
		    debug("Adding Node %d ", neighId);
		  
		
		}// for neighbors

		if(goodNode->parent != -1){
		  next.push_back(goodNode->parent);
		  debug("Also adding this node %d", goodNode->parent );
		}

		if(next.size() == 0){
		  info("we have no more neighbors");
		  cond = false;
		  if(potCm != -1){
		    debug("Last node added belonged to a cm %d, let's see if we should merge them", potCm);
		    (curCand.m_toMerge).push_back(potCm);
		    (tracklets[potCm]->m_toMerge).push_back(curCand.m_id);
		    visitedTracks[curCand.m_id] = 1;
		    visitedTracks[potCm] = 1;
		  }
		} else {

		  debug("Let's add this node and continue");
		  curCand.insertNewNode(goodNode, 1);
		  visited[goodId] = 4;
		  x.erase(x.begin());
		  y.erase(y.begin());
		  x.push_back(goodNode->m_xDet);
		  y.push_back(goodNode->m_yDet);
		}
		
	      }
	      
	      //   cond = false;
	    } // WHILE COND

	  } // right neighbors !
	} // if size if large enough
	
      } // else curcan finished
    } // for tracklet ssize
    }
    info("Number of connected components: %d", tracklets.size());    
  
    for(unsigned int l = 0; l < tracklets.size(); l++){
      PathCandidate &curCand = *(tracklets[l]);
      std::set<int> const *trk = curCand.m_memberIdSet;
      //  if(curCand.m_isValid) {
      std::set<int> *comp = new std::set<int>((*trk));
	    
      connectedComp->push_back(comp);

      //}
    }
    
    int NumConnComp = connectedComp->size();
    ComponentPerEvt.Fill(k, NumConnComp);
    // Store the data for each constructed component
    for(size_t cm = 0 ; cm < connectedComp->size(); ++cm) {
      std::set<int> const* idset = connectedComp->at(cm);
      if(!idset){
	continue;
      }
      std::set<int>::iterator it;
      for( it = idset->begin(); it != idset->end(); ++it) {
     	int detID = *it;// Id of the current detector
	//	printf("CM %d, id %d \n", cm, detID);
     	int d_Index = gr.Find(detID);// Index in the grid
     	GridNode  &node = Ingrid[d_Index];	
	// k = event number, cm = component number,....
	ConnectedCoord.Fill(k, cm, node.m_detID, node.m_x, node.m_y, node.m_z,
			    node.m_xDet, node.m_yDet, node.m_z_Det);
	//printf("%d, %d, %d, %f, %f, %f, %f, %f, %f \n", k, cm, node.m_detID, node.m_x, node.m_y, node.m_z,
	//		    node.m_xDet, node.m_yDet, node.m_z_Det);


      }
    // Print Info
    }

     // ------------------------------------------------------------------------
    #if (WRITE_CM_ASCII_FILE > 0)  
    // Write to an ASCII file
    //////////// Comments to write to file
    std::string header = "cm,x,y,z,mx,my,mz\n";
    //////// End of comments.
    std::ofstream OutTxtFile;
    OutTxtFile.open ("ConnectedComponents.csv");
  
    if (OutTxtFile.is_open()) {
      OutTxtFile << header;
      for(unsigned int l = 0; l < tracklets.size(); l++){
	PathCandidate &curCand = *(tracklets[l]);
	std::vector<int>  *trk = curCand.m_memberList;
	for(size_t cm = 0 ; cm < trk->size(); ++cm) {
	  int detID = trk->at(cm);// Id of the current detector
	  //	printf("CM %d, id %d \n", cm, detID);
	  int d_Index = gr.Find(detID);// Index in the grid
	  GridNode  &node = Ingrid[d_Index];
	  OutTxtFile << l << "," << node.m_x <<"," << node.m_y <<"," << node.m_z <<"," <<	  node.m_xDet<<"," << node.m_yDet<<"," << node.m_z_Det<< '\n';
	  // k = event number, cm = component number,....
	  //	ConnectedCoord.Fill(k, cm, node.m_detID, node.m_x, node.m_y, node.m_z,
	  //			    node.m_xDet, node.m_yDet, node.m_z_Det);
	  //printf("%d, %d, %d, %f, %f, %f, %f, %f, %f \n", k, cm, node.m_detID, node.m_x, node.m_y, node.m_z,
	  //		    node.m_xDet, node.m_yDet, node.m_z_Det);


	}
	// Print Info
      }
      // Writ
    }

  
  
    //OutTxtFile.close();
    #endif

    CollectGridToTree(gr, coord);

    if(connectedComp != 0) {
      for(size_t c = 0; c < connectedComp->size(); ++c) {
        delete connectedComp->at(c);
      }
      delete connectedComp;
    }

  }
  timer.Stop();
  ComponentPerEvt.Write();
  ConnectedCoord.Write();


  
  // Write coordinates ntuple
  coord.Write();
  Out_Put_File.Close();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  std::cout <<"=======================================\n"
            << "Macro finished succesfully.\n"
            << "Real time " << (rtime/totalnumEvt)
	    << " (s/Event), CPU time " << (ctime/totalnumEvt)
            << " (s/Event).\n"
            << '\n';
  }
  // Extract MC track values
  
}
