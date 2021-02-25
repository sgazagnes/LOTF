
#include <iostream>
#include <set>
#include <sstream>
#include <algorithm>
#include  <cmath>
#include <fstream>
#include <vector>
#include <gsl/gsl_poly.h>
#include <cstdlib>
#include <stdlib.h> 
// Root headers
#include "TFile.h"
#include "TNtuple.h"
#include "TStopwatch.h"
#include "TH1.h"

#include "simon_functions.h"
#include "logc.h"
#include "path_queue.h"



bool sortbysec(const pair<int,unsigned short> &a, 
              const pair<int,unsigned short> &b) 
{ 
    return (a.second < b.second); 
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
  dbgfit("Checking curvature = %lf",curv);
  return curv;

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
  float v_y = endTubeB_y - startTubeB_y;
  float v_z = endTubeB_z - startTubeB_z;

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
  dbggrid("Total number of sectors = %d \n Number of left boundary limits = %d\n Number of right boundary limits = %d", SectorLeft.size(), SectorLeft.size(), SectorRight.size());
	   
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
	  currDist =sqrt((CurLftNode.m_x - CurRgtNode.m_x) * (CurLftNode.m_x - CurRgtNode.m_x) +
			 (CurLftNode.m_y - CurRgtNode.m_y) * (CurLftNode.m_y - CurRgtNode.m_y));// +distanceBetweenTube(CurLftNode, CurRgtNode);
	  //if(CurLftNode.m_detID == 1521 || CurRgtNode.m_detID == 1521 )
	  //printf("Distance between node %d and %d is %lf \n",Ingrid[currentLeftSector[l]].m_detID, CurRgtNode.m_detID, currDist);
	  if(currDist <5.)
	    (CurLftNode.m_neighbors).push_back(CurRgtNode.m_detID);

	  
	  /* (CurLftNode.m_x - CurRgtNode.m_x) * (CurLftNode.m_x - CurRgtNode.m_x) +
	    (CurLftNode.m_y - CurRgtNode.m_y) * (CurLftNode.m_y - CurRgtNode.m_y) +
	    (CurLftNode.m_z - CurRgtNode.m_z) * (CurLftNode.m_z - CurRgtNode.m_z);
	  // Update shortest distance

	    if(CurLftNode.m_Layer == CurRgtNode.m_Layer && minDist[1] > currDist) {
	    minDist[1] = currDist;
	    shortestIndex[1] = CurRgtNode.m_detID;//currentRightSector[r];
	  } else if (CurLftNode.m_Layer < CurRgtNode.m_Layer && minDist[0] > currDist){
	    minDist[0] = currDist;
	    shortestIndex[0] = CurRgtNode.m_detID;//currentRightSector[r];
	  } else if (CurLftNode.m_Layer > CurRgtNode.m_Layer && minDist[2] > currDist){
	    minDist[2] = currDist;
	    shortestIndex[2] = CurRgtNode.m_detID;//currentRightSector[r];
	    }*/
	  //limL.Fill(CurLftNode.m_x, CurLftNode.m_y, CurLftNode.m_z);
	  //limR.Fill(CurRgtNode.m_x, CurRgtNode.m_y, CurRgtNode.m_z);
	}// If same layer And not same sector
	

      }// Right limit Loop
      // Found the tube with the shortest distance.
      // if(CurLftNode.m_detID == 1521 || CurLftNode.m_detID == 1521)
      //	printf("%d %d %d  \n",shortestIndex[0],shortestIndex[1],shortestIndex[2]);
      /*	if(shortestIndex[0] > -1 && !(std::find((CurLftNode.m_neighbors).begin(), (CurLftNode.m_neighbors).end(), shortestIndex[0]) != (CurLftNode.m_neighbors).end())){
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

	  }*/
      //  InterSectorPairs.push_back( std::make_pair(currentLeftSector[l], shortestIndex) );
    }// Left limit loop
  }// Sectors loop
  // Now we have all the tube pairs for which we can add virtual nodes
  // between the sectors.
  // DEBUG NTUPLE
  
}





double IntersectionPointSkePar(CoordGrid const &hitMap,
				   GridNode &tubeA, GridNode &tubeB,
				   GridNode &out)
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
  float v_y = endTubeB_y - startTubeB_y;
  float v_z = endTubeB_z - startTubeB_z;

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
    double    a = (double) u_x*u_x + u_y*u_y + u_z*u_z;         // always >= 0
    double    b = (double) u_x*v_x + u_y*v_y + u_z*v_z;
    double    c = (double) v_x*v_x + v_y*v_y + v_z*v_z;         // always >= 0
    double    d = (double) u_x*w_x + u_y*w_y + u_z*w_z; //u.dot(w);
    double    e = (double) v_x*w_x + v_y*w_y + v_z*w_z; // v.dot(w);
    double    D = (double) a*c - b*b;        // always >= 0
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

    double diff_x = (double) ((1.0-sc) * startTubeA_x + sc * endTubeA_x) + ((1.0 - tc) * startTubeB_x + tc * endTubeB_x);
    double diff_y = (double) ((1.0-sc) * startTubeA_y + sc * endTubeA_y) + ((1.0 - tc) * startTubeB_y + tc * endTubeB_y);
    double diff_z = (double) ((1.0-sc) * startTubeA_z + sc * endTubeA_z) + ((1.0 - tc) * startTubeB_z + tc * endTubeB_z);

   double dist = sqrt(pow(((1.0-sc) * startTubeA_x + sc * endTubeA_x) - ((1.0 - tc) * startTubeB_x + tc * endTubeB_x),2)+pow(((1.0-sc) * startTubeA_y + sc * endTubeA_y) - ((1.0 - tc) * startTubeB_y + tc * endTubeB_y),2));
    if(fabs(diff_z/2.) > 140){
      error("Errore, %f", diff_z/2);
      cout << "tube A"<< "\t" << startTubeA_z << "\t" <<   endTubeA_z <<  endl;
      cout << "tube B" << "\t"<< startTubeB_z << "\t" <<  endTubeB_z << endl;
    }
  
    GridNode TransA(tubeA);// Dit kan beter maar voor nu .... FIXME Later


    TransA.m_x     = (tubeA.m_x + tubeB.m_x)/2.0;
    TransA.m_y     = (tubeA.m_y + tubeB.m_y)/2.0;
    TransA.m_xDet  = (tubeA.m_x + tubeB.m_x)/2.0;
    TransA.m_yDet  = (tubeA.m_y + tubeB.m_y)/2.0;
    TransA.m_z     = (float) diff_z/2.;
    TransA.m_z_Det = (float) diff_z/2.;
    // debug("%f", diff_z/2.);
    // Set node as virtual and weight becomes 0 (no added value for
    // the length or area size).
    TransA.m_type        = GridNode::VIRTUAL_NODE;
    TransA.m_weight      = 0;
    TransA.m_SectorLimit = 0;
    TransA.m_LayerLimit  = false;
    TransA.m_halfLength  = 0;// A point. No half length
    // Add parents to the neigboring list
    (TransA.m_neighbors).clear();
    (TransA.m_neighbors).push_back(tubeA.m_detID);
    (TransA.m_neighbors).push_back(tubeB.m_detID);

    
    out = TransA;

    ////////////
    return dist;

}

inline float Det(float a, float b, float c, float d)
{
	return a*d - b*c;
}

bool LineLineIntersect( GridNode &tubeA, GridNode &tubeB, GridNode &tubeC, float &ixOut, float &iyOut, float &izOut) //Output 
{
    //http://mathworld.wolfram.com/Line-LineIntersection.html

  TVector3 dir = tubeC.m_WireDirection;

  float x1 =  tubeA.m_x; 
  float y1 =  tubeA.m_y; 

  float x2 =  tubeB.m_x;
  float y2 =  tubeB.m_y;

  float x3 =  tubeC.m_x - tubeC.m_halfLength * dir[0];
  float y3 =  tubeC.m_y - tubeC.m_halfLength * dir[1];
  float z3 =  tubeC.m_z - tubeC.m_halfLength * dir[2];
  
  float x4 =  tubeC.m_x + tubeC.m_halfLength * dir[0];
  float y4 =  tubeC.m_y + tubeC.m_halfLength * dir[1];
  float z4 =  tubeC.m_z + tubeC.m_halfLength * dir[2];

  float detL1 = Det(x1, y1, x2, y2);
  float detL2 = Det(x3, y3, x4, y4);
  float x1mx2 = x1 - x2;
  float x3mx4 = x3 - x4;
  float y1my2 = y1 - y2;
  float y3my4 = y3 - y4;

  float xnom = Det(detL1, x1mx2, detL2, x3mx4);
  float ynom = Det(detL1, y1my2, detL2, y3my4);
  float denom = Det(x1mx2, y1my2, x3mx4, y3my4);
  
  if(denom == 0.0)//Lines don't seem to cross
    {
      ixOut = NAN;
      iyOut = NAN;
      return false;
      error("Lines do not cross");
    }

  ixOut = xnom / denom;   
  iyOut = ynom / denom;
  izOut = z3 + 2*tubeC.m_halfLength*dir[2]*(iyOut - y3)/(2*tubeC.m_halfLength * dir[1]);
  //if(!isfinite(ixOut) || !isfinite(iyOut)) //Probably a numerical issue
  //    return false;

  return true; //All OK
}

bool PointsLineIntersectLive( GridNode &tubeC, float x1, float x2, float y1, float y2) //Output 
{
    //http://mathworld.wolfram.com/Line-LineIntersection.html

  TVector3 dir = tubeC.m_WireDirection;

  /*float x1 =  tubeA.m_x; 
  float y1 =  tubeA.m_y; 

  float x2 =  tubeB.m_x;
  float y2 =  tubeB.m_y;*/

  float x3 =  tubeC.m_x - tubeC.m_halfLength * dir[0];
  float y3 =  tubeC.m_y - tubeC.m_halfLength * dir[1];
  float z3 =  tubeC.m_z - tubeC.m_halfLength * dir[2];
  
  float x4 =  tubeC.m_x + tubeC.m_halfLength * dir[0];
  float y4 =  tubeC.m_y + tubeC.m_halfLength * dir[1];
  float z4 =  tubeC.m_z + tubeC.m_halfLength * dir[2];

  float detL1 = Det(x1, y1, x2, y2);
  float detL2 = Det(x3, y3, x4, y4);
  float x1mx2 = x1 - x2;
  float x3mx4 = x3 - x4;
  float y1my2 = y1 - y2;
  float y3my4 = y3 - y4;

  float xnom = Det(detL1, x1mx2, detL2, x3mx4);
  float ynom = Det(detL1, y1my2, detL2, y3my4);
  float denom = Det(x1mx2, y1my2, x3mx4, y3my4);

  float distx3, distx4, disx3x4;
  if(denom == 0.0)//Lines don't seem to cross
    {
      //  ixOut = NAN;
      //  iyOut = NAN;
      // error("Lines do not cross,xdet  %f, x3 %f, x4 %f, yDet %f, y3 %f, y4 %f", xnom / denom, x3, x4, ynom / denom, y3, y4);
      return false;
    } else {
    distx3 = sqrt(pow(xnom/denom-x3,2) + pow(ynom/denom -y3,2));
    distx4 = sqrt(pow(xnom/denom-x4,2) + pow(ynom/denom -y4,2));
    disx3x4 = sqrt(pow(x4-x3,2) + pow(y4 -y3,2));
    if(distx3 +distx4 > disx3x4+1){
      // error("Point is not on the line", xnom / denom, x3, x4, ynom / denom, y3, y4);
      return false;
    }
  }
  tubeC.m_xDet = xnom / denom;   
  tubeC.m_yDet = ynom / denom;
  tubeC.m_z_Det = z3 + 2*tubeC.m_halfLength*dir[2]*(tubeC.m_yDet  - y3)/(2*tubeC.m_halfLength * dir[1]);
  //if(!isfinite(ixOut) || !isfinite(iyOut)) //Probably a numerical issue
  //    return false;

  return true; //All OK
}

bool PointsLineIntersectFinal( GridNode &tubeC, float x1, float x2, float y1, float y2) //Output 
{
    //http://mathworld.wolfram.com/Line-LineIntersection.html

  TVector3 dir = tubeC.m_WireDirection;

  /*float x1 =  tubeA.m_x; 
  float y1 =  tubeA.m_y; 

  float x2 =  tubeB.m_x;
  float y2 =  tubeB.m_y;*/

  float x3 =  tubeC.m_x - tubeC.m_halfLength * dir[0];
  float y3 =  tubeC.m_y - tubeC.m_halfLength * dir[1];
  float z3 =  tubeC.m_z - tubeC.m_halfLength * dir[2];
  
  float x4 =  tubeC.m_x + tubeC.m_halfLength * dir[0];
  float y4 =  tubeC.m_y + tubeC.m_halfLength * dir[1];
  float z4 =  tubeC.m_z + tubeC.m_halfLength * dir[2];

  float detL1 = Det(x1, y1, x2, y2);
  float detL2 = Det(x3, y3, x4, y4);
  float x1mx2 = x1 - x2;
  float x3mx4 = x3 - x4;
  float y1my2 = y1 - y2;
  float y3my4 = y3 - y4;

  float xnom = Det(detL1, x1mx2, detL2, x3mx4);
  float ynom = Det(detL1, y1my2, detL2, y3my4);
  float denom = Det(x1mx2, y1my2, x3mx4, y3my4);

  float distx3, distx4, disx3x4;
  if(denom == 0.0)//Lines don't seem to cross
    {
      //  ixOut = NAN;
      //  iyOut = NAN;
      // error("Lines do not cross,xdet  %f, x3 %f, x4 %f, yDet %f, y3 %f, y4 %f", xnom / denom, x3, x4, ynom / denom, y3, y4);
      return false;
    }  else {
    distx3 = sqrt(pow(xnom/denom-x3,2) + pow(ynom/denom -y3,2));
    distx4 = sqrt(pow(xnom/denom-x4,2) + pow(ynom/denom -y4,2));
    disx3x4 = sqrt(pow(x4-x3,2) + pow(y4 -y3,2));
    if(distx3 +distx4 > disx3x4+1){
      if(distx3<distx4){
	tubeC.m_xDet =x3;   
	tubeC.m_yDet = y3;
	tubeC.m_z_Det = z3 + 2*tubeC.m_halfLength*dir[2]*(tubeC.m_yDet  - y3)/(2*tubeC.m_halfLength * dir[1]);
      } else if (distx4<distx3){
	tubeC.m_xDet =x4;   
	tubeC.m_yDet = y4;
	tubeC.m_z_Det = z3 + 2*tubeC.m_halfLength*dir[2]*(tubeC.m_yDet  - y3)/(2*tubeC.m_halfLength * dir[1]);
      }
      // error("Point is not on the line", xnom / denom, x3, x4, ynom / denom, y3, y4);
      return true;
    }
  }
  tubeC.m_xDet = xnom / denom;   
  tubeC.m_yDet = ynom / denom;
  tubeC.m_z_Det = z3 + 2*tubeC.m_halfLength*dir[2]*(tubeC.m_yDet  - y3)/(2*tubeC.m_halfLength * dir[1]);
  //if(!isfinite(ixOut) || !isfinite(iyOut)) //Probably a numerical issue
  //    return false;

  return true; //All OK
}


double IntersectionPointSkeSke(CoordGrid const &hitMap,
				   GridNode &tubeA, GridNode &tubeB,
				   GridNode &out)
{
  TVector3 dirA = tubeA.m_WireDirection;
  float R_A =  tubeA.m_halfLength / sqrt( (dirA[0]*dirA[0]) + (dirA[1]*dirA[1]) + (dirA[2]*dirA[2]) );
  TVector3 dirB = tubeB.m_WireDirection;
  float R_B =  tubeB.m_halfLength / sqrt( (dirB[0]*dirB[0]) + (dirB[1]*dirB[1]) + (dirB[2]*dirB[2]) );

  float startTubeB_x =  tubeB.m_x - tubeB.m_halfLength * dirB[0]-0.10;
  float startTubeB_y =  tubeB.m_y - tubeB.m_halfLength * dirB[1]-0.10;
  float startTubeB_z =  tubeB.m_z - tubeB.m_halfLength * dirB[2];

  float endTubeB_x =  tubeB.m_x + tubeB.m_halfLength * dirB[0]+0.10;
  float endTubeB_y =  tubeB.m_y + tubeB.m_halfLength * dirB[1]+0.10;
  float endTubeB_z =  tubeB.m_z + tubeB.m_halfLength * dirB[2];

  
  float startTubeA_x =  tubeA.m_x - tubeA.m_halfLength * dirA[0]-0.10;
  float startTubeA_y =  tubeA.m_y - tubeA.m_halfLength * dirA[1]-0.10;
  float startTubeA_z =  tubeA.m_z - tubeA.m_halfLength * dirA[2];

  float endTubeA_x =  tubeA.m_x + tubeA.m_halfLength * dirA[0]+0.10;
  float endTubeA_y =  tubeA.m_y + tubeA.m_halfLength * dirA[1]+0.10;
  float endTubeA_z =  tubeA.m_z + tubeA.m_halfLength * dirA[2];
  
  //segment(center - (halflength * direction), center + (halflength * direction))

  float u_x = endTubeA_x - startTubeA_x;
  float u_y = endTubeA_y - startTubeA_y;
  float u_z = endTubeA_z - startTubeA_z;

  float v_x = endTubeB_x - startTubeB_x;
  float v_y = endTubeB_y - startTubeB_y;
  float v_z = endTubeB_z - startTubeB_z;

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
    double    a = (double) u_x*u_x + u_y*u_y + u_z*u_z;         // always >= 0
    double    b = (double) u_x*v_x + u_y*v_y + u_z*v_z;
    double    c = (double) v_x*v_x + v_y*v_y + v_z*v_z;         // always >= 0
    double    d = (double) u_x*w_x + u_y*w_y + u_z*w_z; //u.dot(w);
    double    e = (double) v_x*w_x + v_y*w_y + v_z*w_z; // v.dot(w);
    double    D = (double) a*c - b*b;        // always >= 0
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

    double diff_x = (double) ((1.0-sc) * startTubeA_x + sc * endTubeA_x) + ((1.0 - tc) * startTubeB_x + tc * endTubeB_x);
    double diff_y = (double) ((1.0-sc) * startTubeA_y + sc * endTubeA_y) + ((1.0 - tc) * startTubeB_y + tc * endTubeB_y);
    double diff_z = (double) ((1.0-sc) * startTubeA_z + sc * endTubeA_z) + ((1.0 - tc) * startTubeB_z + tc * endTubeB_z);

    double dist = sqrt(pow(((1.0-sc) * startTubeA_x + sc * endTubeA_x) - ((1.0 - tc) * startTubeB_x + tc * endTubeB_x),2)+pow(((1.0-sc) * startTubeA_y + sc * endTubeA_y) - ((1.0 - tc) * startTubeB_y + tc * endTubeB_y),2));
   
    GridNode TransA(tubeA);// Dit kan beter maar voor nu .... FIXME Later

    // FIXME FIXME Not really optimal (correct??)
     TransA.m_x     = (float)diff_x/2.;
    TransA.m_y     =(float) diff_y/2.;
    TransA.m_xDet  = (float)diff_x/2.;
    TransA.m_yDet  =(float) diff_y/2.;
    TransA.m_z     = (float) diff_z/2.;
    TransA.m_z_Det = (float) diff_z/2.;
    // debug("%f", diff_z/2.);
    // Set node as virtual and weight becomes 0 (no added value for
    // the length or area size).
    TransA.m_type        = GridNode::VIRTUAL_NODE;
    TransA.m_weight      = 0;
    TransA.m_SectorLimit = 0;
    TransA.m_LayerLimit  = false;
    TransA.m_halfLength  = 0;// A point. No half length
    // Add parents to the neigboring list
    (TransA.m_neighbors).clear();

    (TransA.m_neighbors).push_back(tubeA.m_detID);
    (TransA.m_neighbors).push_back(tubeB.m_detID);

    /* tubeA.m_xDet = ((1.0-sc) * startTubeA_x + sc * endTubeA_x);
    tubeA.m_yDet  = (1.0-sc) * startTubeA_y + sc * endTubeA_y;
    tubeA.m_z_Det = (1.0-sc) * startTubeA_z + sc * endTubeA_z;

    tubeB.m_xDet = ((1.0-tc) * startTubeB_x + tc * endTubeB_x);
    tubeB.m_yDet  = (1.0-tc) * startTubeB_y + tc * endTubeB_y;
    tubeB.m_z_Det = (1.0-tc) * startTubeB_z + tc * endTubeB_z;*/

    out = TransA;

    ////////////
    return dist;

}


void IntersectionPointCoord(GridNode &tubeA, GridNode &tubeB)
{
  TVector3 dirA = tubeA.m_WireDirection;
  float R_A =  tubeA.m_halfLength / sqrt( (dirA[0]*dirA[0]) + (dirA[1]*dirA[1]) + (dirA[2]*dirA[2]) );
  TVector3 dirB = tubeB.m_WireDirection;
  float R_B =  tubeB.m_halfLength / sqrt( (dirB[0]*dirB[0]) + (dirB[1]*dirB[1]) + (dirB[2]*dirB[2]) );

  float startTubeB_x =  tubeB.m_x - tubeB.m_halfLength * dirB[0]-0.10;
  float startTubeB_y =  tubeB.m_y - tubeB.m_halfLength * dirB[1]-0.10;
  float startTubeB_z =  tubeB.m_z - tubeB.m_halfLength * dirB[2];

  float endTubeB_x =  tubeB.m_x + tubeB.m_halfLength * dirB[0]+0.10;
  float endTubeB_y =  tubeB.m_y + tubeB.m_halfLength * dirB[1]+0.10;
  float endTubeB_z =  tubeB.m_z + tubeB.m_halfLength * dirB[2];

  
  float startTubeA_x =  tubeA.m_x - tubeA.m_halfLength * dirA[0]-0.10;
  float startTubeA_y =  tubeA.m_y - tubeA.m_halfLength * dirA[1]-0.10;
  float startTubeA_z =  tubeA.m_z - tubeA.m_halfLength * dirA[2];

  float endTubeA_x =  tubeA.m_x + tubeA.m_halfLength * dirA[0]+0.10;
  float endTubeA_y =  tubeA.m_y + tubeA.m_halfLength * dirA[1]+0.10;
  float endTubeA_z =  tubeA.m_z + tubeA.m_halfLength * dirA[2];
  
  //segment(center - (halflength * direction), center + (halflength * direction))

  float u_x = endTubeA_x - startTubeA_x;
  float u_y = endTubeA_y - startTubeA_y;
  float u_z = endTubeA_z - startTubeA_z;

  float v_x = endTubeB_x - startTubeB_x;
  float v_y = endTubeB_y - startTubeB_y;
  float v_z = endTubeB_z - startTubeB_z;

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
  double    a = (double) u_x*u_x + u_y*u_y + u_z*u_z;         // always >= 0
  double    b = (double) u_x*v_x + u_y*v_y + u_z*v_z;
  double    c = (double) v_x*v_x + v_y*v_y + v_z*v_z;         // always >= 0
  double    d = (double) u_x*w_x + u_y*w_y + u_z*w_z; //u.dot(w);
  double    e = (double) v_x*w_x + v_y*w_y + v_z*w_z; // v.dot(w);
  double    D = (double) a*c - b*b;        // always >= 0
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


  //  (TransA.m_neighbors).push_back(tubeA.m_detID);
  //  (TransA.m_neighbors).push_back(tubeB.m_detID);

  tubeA.m_xDet = ((1.0-sc) * startTubeA_x + sc * endTubeA_x);
  tubeA.m_yDet  = (1.0-sc) * startTubeA_y + sc * endTubeA_y;
  tubeA.m_z_Det = (1.0-sc) * startTubeA_z + sc * endTubeA_z;

  tubeB.m_xDet = ((1.0-tc) * startTubeB_x + tc * endTubeB_x);
  tubeB.m_yDet  = (1.0-tc) * startTubeB_y + tc * endTubeB_y;
  tubeB.m_z_Det = (1.0-tc) * startTubeB_z + tc * endTubeB_z;


  ////////////
  return;
}


void Add_VirtualNodes(CoordGrid &hitMap, std::vector < GridNode > &VNodesLayer,  std::vector < GridNode > &VNodesSector)
{
  dbggrid("Computing grid virtual tubes (Neigbor based).");
  // Dummy local counters.
  int NumTubesAdded  = 0;
  // Start node ID for virtual nodes.
  int StartVirtualID = START_VIRTUAL_ID;//6000;

  // List of all nodes available in the detector map.
  std::vector< GridNode > &Ingrid = hitMap.m_grid;
  // Sort all nodes by their layer number
  dbggrid("Sorting nodes, increasing layer number");
  std::sort(Ingrid.begin(), Ingrid.end(), LessThanLayer);
  // Determine the pairs of nodes between the layers for which we want
  // to compute virtual nodes and discard duplicates.
  std::vector< TubeLayerPairProperty > NodePairSet;
  FindNodeBetweenLayerNodePairs(hitMap, NodePairSet);

  // Loop through the possible combinations and create virtual nodes
  // if the actual tubes virtually intersect.
  for(size_t i = 0; i < NodePairSet.size(); ++i) {
    TubeLayerPairProperty const &prop = NodePairSet[i];
    if(prop.isValid) {
      GridNode  &firstNode  = Ingrid[ prop.firstNodeIndex ];
      GridNode  &secondNode = Ingrid[ prop.secondNodeIndex ];
      // Create dummy virtual

      GridNode Dummy_coord;
      // Find intersection point.
      //     if(firstNode.m_type == GridNode::STT_TYPE_SKEW && secondNode.m_type == GridNode::STT_TYPE_SKEW)
      //     //	IntersectionPoint(hitMap, firstNode, secondNode, Dummy_coord) ;
      //     else
      if(firstNode.m_type == GridNode::STT_TYPE_SKEW && secondNode.m_type == GridNode::STT_TYPE_SKEW)
	IntersectionPointSkeSke(hitMap, firstNode, secondNode, Dummy_coord);
      else
	IntersectionPointSkeSke(hitMap, firstNode, secondNode, Dummy_coord);

	// Modify node ID
	//if( firstNode.m_detID == 979 || secondNode.m_detID == 979)
	//  info("%d, %d, %lf", firstNode.m_detID, secondNode.m_detID, distanceBetweenTube(firstNode, secondNode));
	Dummy_coord.m_detID      = StartVirtualID + NumTubesAdded;
	Dummy_coord.m_Orig_detID = Dummy_coord.m_detID;
	std::pair<float, float> r_Theta;
	float theta_deg = Cartesian_To_Polar(Dummy_coord.m_xDet, Dummy_coord.m_yDet, r_Theta);
	Dummy_coord.m_r = r_Theta.first;
	Dummy_coord.m_thetaDeg = theta_deg;

	VNodesLayer.push_back(Dummy_coord);
	NumTubesAdded++;
      // If intersect
    }// END if the pair is valid
  }// End for NodePairSet
  // Reset visiting variable for all nodes in the input graph

  std::vector< std::vector<size_t> > SectorLeft;
  std::vector< std::vector<size_t> > SectorRight;

  for(size_t i = 0; i < 6; ++i) {
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
  dbggrid(" Total number of sectors = %d", SectorLeft.size());
  dbggrid(" Number of left boundary limits = %d" , SectorLeft.size());
  dbggrid(" Number of right boundary limits = %d", SectorRight.size());
  
  double currDist = 0;
  size_t R_Index = 0;
  // Loop all the sectors.
  for(size_t s = 0; s < 6; ++s) { // Numsectors = 6
    std::vector<size_t> const &currentLeftSector = SectorLeft[s];
    // The counter sector is SectorRight[s - 1]. The right limts of
    // the previous sector.
    R_Index = ( s == 0 ) ? (SectorRight.size() - 1) : (s - 1);
    std::vector<size_t> const &currentRightSector = SectorRight[R_Index];
    // For each tube find its nearest neighbor in the other sector.
    for(size_t l = 0; l < currentLeftSector.size(); ++l) {
      GridNode  &CurLftNode = Ingrid[currentLeftSector[l]];
      for(size_t r = 0; r < currentRightSector.size(); ++r) {
	GridNode  &CurRgtNode = Ingrid[currentRightSector[r]];
	// If in the same layer
	if( (labs(CurLftNode.m_Layer - CurRgtNode.m_Layer) <=1 ) &&
	    (CurLftNode.m_Sector != CurRgtNode.m_Sector) &&
	    (CurLftNode.m_type == GridNode::STT_TYPE_SKEW && CurRgtNode.m_type == GridNode::STT_TYPE_SKEW) &&
	    (CurLftNode.m_detID != CurRgtNode.m_detID)
	    ){
	  currDist = distanceBetweenTube(CurLftNode, CurRgtNode);
	  //sqrt((CurLftNode.m_x - CurRgtNode.m_x) * (CurLftNode.m_x - CurRgtNode.m_x) + (CurLftNode.m_y - CurRgtNode.m_y) * (CurLftNode.m_y - CurRgtNode.m_y));
	  
	  //if( CurLftNode.m_detID == 1318 || CurRgtNode.m_detID == 1318)
	  //    /info("%d, %d, %lf", firstNode.m_detID, secondNode.m_detID, distanceBetweenTube(firstNode, secondNode));
	  //info("Node %d and %d, dist %f", CurLftNode.m_detID, CurRgtNode.m_detID, currDist); 
	  if(currDist <3.5){
	    (CurLftNode.m_neighbors).push_back(CurRgtNode.m_detID);
	    /*	    GridNode Dummy_coord;
	    // Find intersection point.
	    IntersectionPointSkeSke(hitMap, CurLftNode, CurRgtNode, Dummy_coord);

	    Dummy_coord.m_detID      = StartVirtualID + NumTubesAdded;
	     Dummy_coord.m_Orig_detID = Dummy_coord.m_detID;
	      std::pair<float, float> r_Theta;
	      float theta_deg = Cartesian_To_Polar(Dummy_coord.m_xDet, Dummy_coord.m_yDet, r_Theta);
	      Dummy_coord.m_r = r_Theta.first;
	      Dummy_coord.m_thetaDeg = theta_deg;

	     VNodesSector.push_back(Dummy_coord);
	     NumTubesAdded++;*/
	  }
	}
      }
    }
  }
  for(size_t i = 0; i < Ingrid.size(); ++i) {
    Ingrid[i].m_visited = false;
  }

  dbggrid(" Determined %d virtual tubes (neighborList)", NumTubesAdded);

}



double IntersectionXY(double startX1, double endX1, double startY1, double endY1, double startX2, double endX2, double startY2, double endY2){
 
  //segment(center - (halflength * direction), center + (halflength * direction))

  double u_x = endX1 - startX1;
  double u_y = endY1 - startY1;

  double v_x = endX2 - startX2;
  double v_y = endY2 - startY2;
  
  double w_x = startX1 - startX2;
  double w_y = startY1 - startY2;

  
  double const SMALL_NUM = std::numeric_limits<double>::epsilon();
  
  double    a = u_x*u_x + u_y*u_y;
  double    b = u_x*v_x + u_y*v_y;
  double    c = v_x*v_x + v_y*v_y;
  double    d = u_x*w_x + u_y*w_y;
  double    e = v_x*w_x + v_y*w_y;
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
    
    sc = (std::abs(sN) < SMALL_NUM ? 0.0 : sN / sD);

    return sc;

}
