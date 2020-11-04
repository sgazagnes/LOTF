
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

#include "simon_functions.h"

int determineSkewed_XYPlane_new( CoordGrid &hitMap, GridNode const &VNode,
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
      if(  SK_NeighNode.m_type == GridNode::STT_TYPE_SKEW ) {
        // Current or next layer
        if( SK_NeighNode.m_Layer == LocalCurrentLayer ||
            SK_NeighNode.m_Layer == LocalNextLayer ) {
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
      else if( SK_NeighNode.m_type == GridNode::VIRTUAL_NODE ) {
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


/* sortNeighbors */

bool sortNeighbors(CoordGrid &gr, GridNode *currentNode, std::vector<int> &prev, std::vector<int> &same, std::vector<int> &next, std::vector<int> &virt, char *visited, int *dir){

  int curDir = *dir;
  std::vector< GridNode > &Ingrid  = gr.m_grid;
  int curLayer 	= currentNode->m_Layer;
  int curId   	= currentNode->m_detID;
  bool cond 	= true;

  
  for(int i = 0; i < currentNode->m_neighbors.size(); i++){
    int neighId 	= currentNode->m_neighbors[i];
    int neighIdx 	= gr.Find(neighId);
    GridNode &neighNode = Ingrid[neighIdx];
    
    // debug("Node %d has one neig %d", curId, neighId);
    
    if(neighNode.m_type == GridNode::VIRTUAL_NODE){
      virt.push_back(neighId);
      // if(neighId == 12454)
//	debug("%d %d", neighNode.m_neighbors[0],  neighNode.m_neighbors[1]);
      //   determineSkewed_XYPlane_perso( *gr, *neighNode, ListOfSkewedNodesIndex, ListOfVirtualNodesIndex,visited);
      //error("%d curlaye %d, virtual layer %d", neighId, curLayer, neighNode->m_Layer);
      continue;
    }	 

    if(visited[neighId] == 0){ // Not visited
      if(neighNode.m_Layer > curLayer){
	//	debug("Node %d has one neigh up %d", curId, neighId);
	next.push_back(neighId);
	curDir |= UP;
	visited[neighId] = 2;
      }

      else if( neighNode.m_Layer < curLayer) {
	//	debug("Node %d has one neigh down %d", curId, neighId);
	prev.push_back(neighId);
	curDir |= DOWN;
	visited[neighId] = 2;
      }

      else {
	//	debug("Node %d has one neigh on the same %d", curId, neighId);
	same.push_back(neighId);
	visited[neighId] = 2;
      }
    } // if not visited

    else if(visited[neighId] == 4){ // Already in CM
      //   debug("Node %d has already been connected, tricky", neighId);
      cond = false;
      if(!(std::find(same.begin(), same.end(), neighId) != same.end()))
	same.push_back(neighId);

    }

    //  else
      // debug("Node %d has already been added to the list", neighId);
  }

  	
  if( curDir > 4 || same.size() > 1){
    info("A complex case to be solved later");
    cond = false;
  }
  
  *dir = curDir;
  //debug("Neig cond %d", cond);
  return cond;
}


/* resetLists */


void resetLists(char *visited, std::vector<int> &prev, std::vector<int> &same, std::vector<int> &next){
  for (int i = 0; i < same.size(); i++)
    if(visited[same[i]] < 4)
      visited[same[i]] =0;
	      
  for (int i = 0; i < next.size(); i++)
    if(visited[next[i]] < 4)
      visited[next[i]] =0;

  for (int i = 0; i < prev.size(); i++)
    if(visited[prev[i]] < 4)
      visited[prev[i]] =0;
}


/* removeIdFromNeigh */


void removeIdFromNeigh(GridNode *neighNode, std::vector<int> *prevNodes, int curId){
  for(int i = 0; i < prevNodes->size(); i++){
    if(prevNodes->at(i) != neighNode->m_detID){ 		
      (neighNode->m_neighbors).erase(std::remove((neighNode->m_neighbors).begin(), (neighNode->m_neighbors).end(),prevNodes->at(i)), (neighNode->m_neighbors).end());
      if(i == 0)
	neighNode->parent = curId;
    }
  }
}


/* areAdjacent */


bool areAdjacent(CoordGrid &gr, std::vector<int> *v){
  int adjacent = 0;
  std::vector< GridNode > &Ingrid = gr.m_grid;

  for (int i = 0; i < v->size(); i++){
    int neighId   = v->at(i);
    int neighIdx  = gr.Find(neighId);
    GridNode &neighNode = Ingrid[neighIdx];
    //prevNodes.push_back(neighId);
    
    for (int j = i+1; j < v->size (); j++){
      //   debug("Are %d and %d connected?", neighId, v->at(j));
      if(neighId == v->at(j)) error("areAdjacent: should not be in vector");
      else if(neighNode.IsNeighboring(v->at(j)))
	adjacent++;
      // else
    }
  }
  
  if(adjacent >= v->size() -1)
    return true;
  else
    return false;
}


/* polyFit */

double *polyFit(std::vector<double>  x, std::vector<double>  y, int n){
  int i,j,k,N;
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
   
  return a;
}

/* returnDirection */ 

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

double returnAngle(double x1, double x2, double x3, double y1, double y2, double y3){
  
 
  float abx =  x2 - x1;
  float aby =  y2 - y1 ;
  float cbx =  x2 - x3;
  float cby =  y2 - y3 ;

  float dot = (abx * cbx + aby * cby); // dot product
  float cross = (abx * cby - aby * cbx); // cross product

  float alpha = atan2(cross, dot);
  float angleDeg = alpha * 180. / 3.14159265 ;
  // return (int) floor(alpha * 180. / pi + 0.5);
  return angleDeg;

}

double returnCurvature(double x1, double x2, double x3, double y1, double y2, double y3){
  
  double len1 = sqrt(pow(x1-x3,2)+pow(y1-y3,2));
  double len2 = sqrt(pow(x2-x3,2)+pow(y2-y3,2));
  double len3 = sqrt(pow(x1-x2,2)+pow(y1-y2,2));
  double area = fabs(x1*(y2 - y3)+ x2*(y3-y1) + x2*(y1-y2))/2.;
  double curv = area/(len1*len2*len3);
  debug("Checking curvature = %lf",curv);
  return curv;

}



double nodeDistanceToLinearFit(double xdet, double ydet, double *x_coef, double *y_coef){
  
  //double xdet = (double) node->m_r/  sqrt( 2*pow(40,2));//node->m_xDet;
  //double ydet = (double)  (node->m_thetaDeg+180.) /360.;

  double newx_coef[2] = {x_coef[0] - xdet, x_coef[1]};
  double newy_coef[2] = {y_coef[0] - ydet, y_coef[1]};

  double vectortanx[2] = {x_coef[1], 0.};
  double vectortany[2] = {y_coef[1], 0.};
	      
  double d[3] = {newx_coef[0]*vectortanx[0] + newy_coef[0]*vectortany[0],
		 newx_coef[0]*vectortanx[1] + newx_coef[1]*vectortanx[0] +
		 newy_coef[0]*vectortany[1] + newy_coef[1]*vectortany[0],
		 newx_coef[1]*vectortanx[1] + newy_coef[1]*vectortany[1]  };

  double x0[2];
  int nroot = gsl_poly_solve_quadratic(d[2], d[1], d[0], x0, x0+1);

  double xIntersect, yIntersect, currDist;
      
  if(nroot ==1){
    xIntersect = gsl_poly_eval(x_coef, 2, x0[0]);
    yIntersect = gsl_poly_eval(y_coef, 2, x0[0]);
    currDist = sqrt(pow(xIntersect - xdet,2) + pow(yIntersect - ydet,2));
  }

		
  for (int j = 0; j < nroot; j++){
    double newx = gsl_poly_eval(x_coef, 2, x0[j]);
    double newy = gsl_poly_eval(y_coef, 2, x0[j]);
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

  return currDist;
}


double nodeDistanceToQuadFit(double xdet, double ydet, double *x_coef, double *y_coef){

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
  return currDist;
}


/* fitNextId */

int fitNextId(CoordGrid &gr, PathCandidate &cand, std::vector<int> &next, int k){
  
  std::vector< GridNode > &Ingrid = gr.m_grid;
  int goodId     = -1;
  int method;
  int degree, nElts;
  std::vector<double> x =  std::vector<double>( cand.m_x );
  std::vector<double> y =  std::vector<double>( cand.m_y ); 
  std::vector<double> r =  std::vector<double>( cand.m_r );
  std::vector<double> theta =  std::vector<double>( cand.m_theta ); 

  for (int i = 0; i <x.size(); i++){
    r[i] /= sqrt( 2*pow(40,2));
    theta[i] = (theta[i]+180.) /360.;    
  }
  double prevtheta = theta.back();

  //x.erase(std::remove(x.begin(), x.end(), -1), x.end());
  // y.erase(std::remove(y.begin(), y.end(), -1), y.end());
  // r.erase(std::remove(r.begin(), r.end(), -1), r.end());
  // theta.erase(std::remove(theta.begin(), theta.end(), -1), theta.end());
  
  if(k == 0){
    std::reverse(x.begin(),x.end());
    std::reverse(y.begin(),y.end());
    std::reverse(r.begin(),r.end());
    std::reverse(theta.begin(),theta.end());
  }

  //debug("Fast-checking unphysical neighbors");

  std::vector<int> plausible;
  std::vector<int> uncertain;
  std::vector<int> unlikely;
  std::vector<int> *tocheck;

  
  float tol = 90.;
  for (int i = 0; i <next.size(); i++){
    int curId = next[i];
    int curIdx = gr.Find(curId);
    GridNode *node = &Ingrid[curIdx];
    double xdet = node->m_x;
    double ydet = node->m_y;
    double rdet = (double) node->m_r /  sqrt( 2*pow(40,2));
    double thetadet = (double) (node->m_thetaDeg+180.) /360.;
    if(node->m_type == GridNode::STT_TYPE_SKEW)
      tol = 50;
    if(thetadet > 0.85 && prevtheta < 0.15)
      thetadet -= 1;
    else if(thetadet < 0.15 && prevtheta > 0.85)
      thetadet += 1.;
    //  debug("Points %lf, %lf \t %lf, %lf \t %lf, %lf",r[r.size()-2], r[r.size()-1], rdet, theta[theta.size()-2], theta[theta.size()-1], thetadet);
    float angle_r = returnAngle(r[r.size()-2], r[r.size()-1], rdet, theta[theta.size()-2], theta[theta.size()-1], thetadet);
    
    float angle_xy = returnAngle(x[x.size()-2], x[x.size()-1], xdet, y[y.size()-2], y[y.size()-1], ydet);
    debug("Angle with %d is %f, %f", curId, angle_r, angle_xy);
      if(fabs(angle_r) > 90 && fabs(angle_xy) > tol)
      plausible.push_back(curId);
      else if((fabs(angle_r) > 90 || fabs(angle_xy) > tol) && fabs(angle_r) > 50 && fabs(angle_xy) > 50)
	uncertain.push_back(curId);
      else
	unlikely.push_back(curId);
  }

  if(plausible.size() == 1){
    goodId = plausible[0];
    info("Only one good choice %d", goodId);
    return goodId;
  }else  if(plausible.size() > 0){
    info("We found %d promising candidates", plausible.size());
    tocheck = &plausible;
  }  else if (uncertain.size() > 0) {
    if(uncertain.size() == 1){
      
      goodId = uncertain[0];
      info("Only one good choice %d", goodId);
      return goodId;
    }
    info("No promising, but still possible with %d cand", uncertain.size());
    //return -1;
    tocheck = &uncertain;
  } else {
    info("Unless we want to go backwards, we should stop");
    return -1;
  }

  // Checking Track angle in polar coord 

  // debug("Points %lf, %lf \t %lf, %lf \t %lf, %lf",r[0], theta[0], r[r.size()/2], theta[theta.size()/2],r[r.size()-1], theta[theta.size()-1]);
  
  // double curv = returnCurvature(x[0], x[x.size()/2], x[x.size()-1], y[0], y[y.size()/2], y[y.size()-1]);

  float angle = returnAngle(r[0], r[r.size()/2], r[r.size()-1],  theta[0], theta[theta.size()/2], theta[theta.size()-1]);
  //debug("Angle of track is %f", angle);
  

  if(fabs(angle) < 170){
    info("Quadratic Fit");
    method = 1;
    degree = 2;
    nElts  = x.size() -1;
  } else {
    info("Linear Fit");
    method = 0;
    degree = 1;
    nElts = x.size()-1;
  }
 
  std::vector<double> p;
  double minDist = std::numeric_limits<double>::max();
    
  p.push_back(0.);
  for (int i = 0; i <  nElts; i++){ 
    double newval = p[i] + sqrt(pow(r[i+1]-r[i],2.) + pow(theta[i+1]-theta[i],2.));
    p.push_back(newval);
  }

  double *x_coef = polyFit(p, r, degree);
  double *y_coef = polyFit(p, theta, degree);

  
  //  if(method == 0){ // linear

    
  
  for (int i = 0; i < tocheck->size(); i++){
    int curId = tocheck->at(i);
    int curIdx = gr.Find(curId);
    GridNode *node = &Ingrid[curIdx];
    double xdet = (double) node->m_r/  sqrt( 2*pow(40,2));//node->m_xDet;
    double ydet = (double)  (node->m_thetaDeg+180.) /360.;     
    if(ydet > 0.85 && prevtheta < 0.15)
      ydet -= 1;
    else if(ydet < 0.15 && prevtheta > 0.85)
      ydet += 1.;
    double currDist = method == 0? nodeDistanceToLinearFit(xdet, ydet, x_coef, y_coef): nodeDistanceToQuadFit(xdet, ydet, x_coef, y_coef);      

     debug("Node %d is at %lf", curId, currDist);
		
    if(minDist > currDist){
      minDist = currDist;		  
      goodId = curId;	
    }        
  }
  
  return goodId;

  if(goodId != -1) info("The good id is %d, and is at distance %lf", goodId,minDist);
  else info("No good ID found");
  

}


/*distanceBetweenTuves */

double distanceBetweenTube(GridNode & tubeA, GridNode & tubeB) 
{

  TVector3 dirA = tubeA.m_WireDirection;
  float R_A =  tubeA.m_halfLength / sqrt( (dirA[0]*dirA[0]) + (dirA[1]*dirA[1]) + (dirA[2]*dirA[2]) );
  TVector3 dirB = tubeB.m_WireDirection;
  float R_B =  tubeB.m_halfLength / sqrt( (dirB[0]*dirB[0]) + (dirB[1]*dirB[1]) + (dirB[2]*dirB[2]) );

  float startTubeB_x =  tubeB.m_x - tubeB.m_halfLength * dirB[0];
  float startTubeB_y =  tubeB.m_y - tubeB.m_halfLength * dirB[1];
  float startTubeB_z =  tubeB.m_z - tubeB.m_halfLength * dirB[2];

  float endTubeB_x =  tubeB.m_x + tubeB.m_halfLength * dirB[0];
  float endTubeB_y =  tubeB.m_y + tubeB.m_halfLength * dirB[1];
  float endTubeB_z =  tubeB.m_z + tubeB.m_halfLength * dirB[2];

  
  float startTubeA_x =  tubeA.m_x - tubeA.m_halfLength * dirA[0];
  float startTubeA_y =  tubeA.m_y - tubeA.m_halfLength * dirA[1];
  float startTubeA_z =  tubeA.m_z - tubeA.m_halfLength * dirA[2];

  float endTubeA_x =  tubeA.m_x + tubeA.m_halfLength * dirA[0];
  float endTubeA_y =  tubeA.m_y + tubeA.m_halfLength * dirA[1];
  float endTubeA_z =  tubeA.m_z + tubeA.m_halfLength * dirA[2];
  
  //segment(center - (halflength * direction), center + (halflength * direction))

  float u_x = endTubeA_x - startTubeA_x;
  float u_y = endTubeA_y - startTubeA_y;
  float u_z = endTubeA_z - startTubeA_z;

  float v_x = endTubeB_x - startTubeB_x;
  float v_y = endTubeA_y - startTubeA_y;
  float v_z = endTubeA_z - startTubeA_z;

  float w_x = startTubeA_x - startTubeB_x;
  float w_y = startTubeA_y - startTubeB_y;
  float w_z = startTubeA_z - startTubeB_z;
  
  /* GRVector3 P0 = start;
    GRVector3 P1 = end;
    GRVector3 Q0 = line.start;
    GRVector3 Q1 = line.end;*/

    double const SMALL_NUM = std::numeric_limits<double>::epsilon();
    /*   GRVector3   u = P1 - P0;
    GRVector3   v = Q1 - Q0;
    GRVector3   w = P0 - Q0;*/
    double    a = u_x*u_x + u_y*u_y + u_z*u_z;         // always >= 0
    double    b = u_x*v_x + u_y*v_y + u_z*v_z;
    double    c = v_x*v_x + v_y*v_y + v_z*v_z;         // always >= 0
    double    d = u_x*w_x + u_y*w_y + u_z*w_z; //u.dot(w);
    double    e = v_x*w_x + v_y*w_y + v_z*w_z; // v.dot(w);
    double    D = a*c - b*b;        // always >= 0
    double    sc, sN, sD = D;       // sc = sN / sD, default sD = D >= 0
    double    tc, tN, tD = D;       // tc = tN / tD, default tD = D >= 0

    // compute the line parameters of the two closest points
    if (D < SMALL_NUM) { // the lines are almost parallel
        sN = 0.0;         // force using point P0 on segment S1
        sD = 1.0;         // to prevent possible division by 0.0 later
        tN = e;
        tD = c;
    }
    else {                 // get the closest points on the infinite lines
        sN = (b*e - c*d);
        tN = (a*e - b*d);
        if (sN < 0.0) {        // sc < 0 => the s=0 edge is visible
            sN = 0.0;
            tN = e;
            tD = c;
        }
        else if (sN > sD) {  // sc > 1  => the s=1 edge is visible
            sN = sD;
            tN = e + b;
            tD = c;
        }
    }

    if (tN < 0.0) {            // tc < 0 => the t=0 edge is visible
        tN = 0.0;
        // recompute sc for this edge
        if (-d < 0.0)
            sN = 0.0;
        else if (-d > a)
            sN = sD;
        else {
            sN = -d;
            sD = a;
        }
    }
    else if (tN > tD) {      // tc > 1  => the t=1 edge is visible
        tN = tD;
        // recompute sc for this edge
        if ((-d + b) < 0.0)
            sN = 0;
        else if ((-d + b) > a)
            sN = sD;
        else {
            sN = (-d + b);
            sD = a;
        }
    }
    // finally do the division to get sc and tc
    sc = (std::abs(sN) < SMALL_NUM ? 0.0 : sN / sD);
    tc = (std::abs(tN) < SMALL_NUM ? 0.0 : tN / tD);

    //  GRVector3 diff = ((1.0 - sc) * P0 + sc * P1) - ((1.0 - tc) * Q0 + tc * Q1);
    double diff_x = ((1.0-sc) * startTubeA_x + sc * endTubeA_x) - ((1.0 - tc) * startTubeB_x + tc * endTubeB_x);
    double diff_y = ((1.0-sc) * startTubeA_y + sc * endTubeA_y) - ((1.0 - tc) * startTubeB_y + tc * endTubeB_y);
    double diff_z = ((1.0-sc) * startTubeA_z + sc * endTubeA_z) - ((1.0 - tc) * startTubeB_z + tc * endTubeB_z);

    return sqrt(diff_x*diff_x + diff_y*diff_y + diff_z*diff_z);
}


/* Fix_InterSector_Nodes */

void Fix_InterSector_Nodes(CoordGrid &hitMap, size_t const numSectors)
{
  // Fetch the list of all STT tubes.
  std::vector< GridNode > &Ingrid = hitMap.m_grid;
 
  // Sector Counting starts from zero.
  // List of per sector tube indices in the map
  std::vector< std::vector<size_t> > SectorLeft;
  std::vector< std::vector<size_t> > SectorRight;

  for(size_t i = 0; i < numSectors; ++i) {
    SectorLeft.push_back(std::vector<size_t>());
    SectorRight.push_back(std::vector<size_t>());
  }
  // DEBUG
  TNtuple limL ("SelectedSectorLeft","SELECTED SKEWED BETWEEN SECTORS","x:y:z");
  TNtuple limR ("SelectedSectorRight","SELECTED SKEWED BETWEEN SECTORS","x:y:z");
  // DEBUG
  // Find sector boundaries per sector and per side
  for(size_t i = 0; i < Ingrid.size(); ++i) {
    GridNode const &tube = Ingrid[i];
    if( (tube.m_type != GridNode::VIRTUAL_NODE) ) {
      // Left or right boundary limit
      if(tube.m_SectorLimit == -1) {
	(SectorLeft[tube.m_Sector]).push_back(i);
	limL.Fill(tube.m_x, tube.m_y, tube.m_z);
      }
      // Right limit
      else if(tube.m_SectorLimit == 1) {
	(SectorRight[tube.m_Sector]).push_back(i);
	limR.Fill(tube.m_x, tube.m_y, tube.m_z);
      }
    }
  }
  //_______
  std::cout << "<DEBUG> Total number of sectors = " << SectorLeft.size()
	    << "\n\t<-I-> Number of left boundary limits = " << SectorLeft.size()
	    << "\n\t<-I-> Number of right boundary limits = " << SectorRight.size() 
	    << '\n';
  /*
   *   0 / \ 5
   *   1 | | 4
   *   2 \ / 3
   */
  double minDist[3];// = std::numeric_limits<double>::max();
  double currDist = 0;
  // Index of the tube with the shortest dstance
  int shortestIndex[3];
  size_t R_Index = 0;
  // Loop all the sectors.
  for(size_t s = 0; s < numSectors; ++s) {
    std::vector<size_t> const &currentLeftSector = SectorLeft[s];
    // The counter sector is SectorRight[s - 1]. The right limts of
    // the previous sector.
    R_Index = ( s == 0 ) ? (SectorRight.size() - 1) : (s - 1);
    std::vector<size_t> const &currentRightSector = SectorRight[R_Index];
    // For each tube find its nearest neighbor in the other sector.
    for(size_t l = 0; l < currentLeftSector.size(); ++l) {
      GridNode  &CurLftNode = Ingrid[currentLeftSector[l]];
      minDist[0] = std::numeric_limits<double>::max();
      minDist[1] = std::numeric_limits<double>::max();
      minDist[2] = std::numeric_limits<double>::max();
      shortestIndex[0] = -1;
      shortestIndex[1] = -1;
      shortestIndex[2] = -1;
      //  info("Node %d", CurLftNode.m_detID);
      for(size_t r = 0; r < currentRightSector.size(); ++r) {
	GridNode  &CurRgtNode = Ingrid[currentRightSector[r]];
	// If in the same layer
	if( (labs(CurLftNode.m_Layer - CurRgtNode.m_Layer) <=1 ) &&
	    (CurLftNode.m_Sector != CurRgtNode.m_Sector) &&
	    (CurLftNode.m_type == GridNode::STT_TYPE_SKEW && CurRgtNode.m_type == GridNode::STT_TYPE_SKEW) &&
	    (CurLftNode.m_detID != CurRgtNode.m_detID)
	  ){
	  currDist =distanceBetweenTube(CurLftNode, CurRgtNode);
	  /*(CurLftNode.m_x - CurRgtNode.m_x) * (CurLftNode.m_x - CurRgtNode.m_x) +
	    (CurLftNode.m_y - CurRgtNode.m_y) * (CurLftNode.m_y - CurRgtNode.m_y) +
	    (CurLftNode.m_z - CurRgtNode.m_z) * (CurLftNode.m_z - CurRgtNode.m_z);*/
	  // Update shortest distance
	  // if(Ingrid[currentLeftSector[l]].m_detID == 1388 ||CurRgtNode.m_detID == 1560 )
	  //   printf("Distance between node %d and %d is %lf \n",Ingrid[currentLeftSector[l]].m_detID, CurRgtNode.m_detID, currDist);
	  if(CurLftNode.m_Layer == CurRgtNode.m_Layer && minDist[1] > currDist) {
	    minDist[1] = currDist;
	    shortestIndex[1] = CurRgtNode.m_detID;//currentRightSector[r];
	  } else if (CurLftNode.m_Layer < CurRgtNode.m_Layer && minDist[0] > currDist){
	    minDist[0] = currDist;
	    shortestIndex[0] = CurRgtNode.m_detID;//currentRightSector[r];
	  } else if (CurLftNode.m_Layer > CurRgtNode.m_Layer && minDist[2] > currDist){
	    minDist[2] = currDist;
	    shortestIndex[2] = CurRgtNode.m_detID;//currentRightSector[r];
	  }
	  //limL.Fill(CurLftNode.m_x, CurLftNode.m_y, CurLftNode.m_z);
	  //limR.Fill(CurRgtNode.m_x, CurRgtNode.m_y, CurRgtNode.m_z);
	}// If same layer And not same sector
	

      }// Right limit Loop
      // Found the tube with the shortest distance.
	if(shortestIndex[0] > -1 && !(std::find((CurLftNode.m_neighbors).begin(), (CurLftNode.m_neighbors).end(), shortestIndex[0]) != (CurLftNode.m_neighbors).end())){
	  (CurLftNode.m_neighbors).push_back(shortestIndex[0]);
	  // debug("New neighbors %d and %d, %lf", CurLftNode.m_detID, shortestIndex[0], minDist[0]);
	}
	if(shortestIndex[1] > -1 && !(std::find((CurLftNode.m_neighbors).begin(), (CurLftNode.m_neighbors).end(), shortestIndex[1]) != (CurLftNode.m_neighbors).end())){
	  (CurLftNode.m_neighbors).push_back(shortestIndex[1]);
	  //  debug("New neighbors %d and %d, %lf", CurLftNode.m_detID, shortestIndex[1], minDist[1]);
	}
	if(shortestIndex[2] > -1 && !(std::find((CurLftNode.m_neighbors).begin(), (CurLftNode.m_neighbors).end(), shortestIndex[2]) != (CurLftNode.m_neighbors).end())){
	  (CurLftNode.m_neighbors).push_back(shortestIndex[2]);
	  //  debug("New neighbors %d and %d, %lf", CurLftNode.m_detID, shortestIndex[2], minDist[2]);

	}
      //  InterSectorPairs.push_back( std::make_pair(currentLeftSector[l], shortestIndex) );
    }// Left limit loop
  }// Sectors loop
  // Now we have all the tube pairs for which we can add virtual nodes
  // between the sectors.
  // DEBUG NTUPLE
  
}


void addTracklets (CoordGrid &gr, PathCandidate *newCand, PathCandidate &mergeCand,  int curdir, int mergedir){

  std::vector< GridNode > &Ingrid = gr.m_grid;

  
  if(mergedir == 1){
    for(int i =  (mergeCand.m_memberList)->size()-1; i >=0; i--){
      int curid = (mergeCand.m_memberList)->at(i);
      //   debug("adding id %d", curid);
      if(newCand->isInCandidate(curid)) continue;
      int curidx = gr.Find(curid);
      GridNode* node = &Ingrid[curidx];
      newCand->insertNewNode(gr,node, curdir == 1? newCand->m_memberList->end(): newCand->m_memberList->begin());
    }
    
  } else {

    for(int i =  0; i < mergeCand.m_memberList->size(); i++){
      int curid = mergeCand.m_memberList->at(i);
      // debug("adding id %d", curid);

      if(newCand->isInCandidate(curid)) continue;
      int curidx = gr.Find(curid);
      GridNode* node = &Ingrid[curidx];
      newCand->insertNewNode(gr, node, curdir == 1? newCand->m_memberList->end(): newCand->m_memberList->begin());
    }
  }
  
  mergeCand.m_isMerged = 1;
  mergeCand.m_isValid = 0;

}
