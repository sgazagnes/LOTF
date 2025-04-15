
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
//#include <opencv2/opencv.hpp>

// Root headers
#include "TFile.h"
#include "TNtuple.h"
#include "TStopwatch.h"
#include "TH1.h"

#include "auxfunctions.h"
#include "logc.h"
//#include "path_queue.h"


#define _USE_MATH_DEFINES
#include <math.h>

#include <random>
#include "circle.h"


bool sortbysec(const pair<int,unsigned short> &a, 
              const pair<int,unsigned short> &b) 
{ 
    return (a.second < b.second); 
} 


inline float Det(float a, float b, float c, float d)
{
	return a*d - b*c;
}


double Cartesian_To_Polar(float const x, float const y,
                          std::pair<float,float>& polarOut,
                          bool useSign)
{
  float radius = sqrt( x*x + y*y);
  
  float theta = 0;
  
  /* Expressed in radians. Taken into account the sign of the
     coordinates.*/
  if(useSign) {
    theta = atan2(y, x);
  }
  else {// Do not take into account the sign
    // x == 0
    if( !(x > 0.0) && !(x < 0.0) ){
      theta = (M_PI / 2.0);
    }
    else if((y/x) >= 0 ) {
      theta = atan(y/x);
    }
    else {
      theta = M_PI + (atan(y/x));
    }
  }
  
  polarOut.first  = radius;
  polarOut.second = theta;
  
  // Return theta in degrees.
  return ( (theta * 180) / M_PI);// 
}


//______________________________ GridToNtuple    ____________________________
TNtuple* GridToNtuple(std::vector < GridNode > const &VNodes, std::string const &name)
{
  TNtuple* out = new TNtuple(name.c_str(),"Grid To Ntuple","x:y:det_z:z");
  for(size_t i = 0; i < VNodes.size(); ++i) {
    GridNode const &tube = VNodes[i];
    out->Fill(tube.m_x, tube.m_y, tube.m_zDet, tube.m_z);
  }
  return out;
}
//_________________________ END GridToNtuple ________________________________

/* polyFit */

double *polyFit(std::vector<double>  x, std::vector<double>  y, int n){
  int i,j,k,N;
  double X[2*n+1];        //Array that will store the values of sigma(xi),sigma(xi^2),sigma(xi^3)....sigma(xi^2n)
  for (i=0;i<2*n+1;i++)
    {
      X[i]=0;
      for (j=0;j<(int) x.size();j++)
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
      for (j=0;j<(int) x.size();j++)
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
  float aby =  y2 - y1;
  float cbx =  x2 - x3;
  float cby =  y2 - y3;

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

double distanceBetweenTube(GridNode &tubeA, GridNode &tubeB) 
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




//IntersectionPointSkePar: [x,y] = Middle, z = intersect
double IntersectionPointTubes(GridNode &tubeA, GridNode &tubeB, GridNode &out, int xymid)
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
  
  float u_x = endTubeA_x - startTubeA_x;
  float u_y = endTubeA_y - startTubeA_y;
  float u_z = endTubeA_z - startTubeA_z;

  float v_x = endTubeB_x - startTubeB_x;
  float v_y = endTubeB_y - startTubeB_y;
  float v_z = endTubeB_z - startTubeB_z;

  float w_x = startTubeA_x - startTubeB_x;
  float w_y = startTubeA_y - startTubeB_y;
  float w_z = startTubeA_z - startTubeB_z;
  

  double const SMALL_NUM = std::numeric_limits<double>::epsilon();
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

  double diff_x = (double) ((1.0-sc) * startTubeA_x + sc * endTubeA_x)
    + ((1.0 - tc) * startTubeB_x + tc * endTubeB_x);
  double diff_y = (double) ((1.0-sc) * startTubeA_y + sc * endTubeA_y)
    + ((1.0 - tc) * startTubeB_y + tc * endTubeB_y);
  double diff_z = (double) ((1.0-sc) * startTubeA_z + sc * endTubeA_z)
    + ((1.0 - tc) * startTubeB_z + tc * endTubeB_z);
  double dist = sqrt(pow(((1.0-sc) * startTubeA_x + sc * endTubeA_x) - ((1.0 - tc) * startTubeB_x + tc * endTubeB_x),2)+pow(((1.0-sc) * startTubeA_y + sc * endTubeA_y) - ((1.0 - tc) * startTubeB_y + tc * endTubeB_y),2));
  
  if(fabs(diff_z/2.) > 140){
    error("Errore, %f", diff_z/2);
    cout << "tube A"<< "\t" << startTubeA_z << "\t" <<  endTubeA_z << endl;
    cout << "tube B" << "\t"<< startTubeB_z << "\t" <<  endTubeB_z << endl;
  }
  
  GridNode TransA(tubeA);
    
  TransA.m_x     = xymid ? (tubeA.m_x + tubeB.m_x)/2.0: (float) diff_x/2.;
  TransA.m_y     = xymid ? (tubeA.m_y + tubeB.m_y)/2.0: (float) diff_y/2.;
  TransA.m_xDet  = xymid ? (tubeA.m_x + tubeB.m_x)/2.0: (float) diff_x/2.;
  TransA.m_yDet  = xymid ? (tubeA.m_y + tubeB.m_y)/2.0: (float) diff_y/2.;
  TransA.m_z     = (float) diff_z/2.;
  TransA.m_zDet = (float) diff_z/2.;

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

// Intersection line with no correction of coordinates
bool PointsLineIntersectLive( GridNode &tubeC, float x1, float x2, float y1, float y2) //Output 
{
    //http://mathworld.wolfram.com/Line-LineIntersection.html

  TVector3 dir = tubeC.m_WireDirection;

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
  tubeC.m_zDet = z3 + 2*tubeC.m_halfLength*dir[2]*(tubeC.m_yDet  - y3)/(2*tubeC.m_halfLength * dir[1]);
  //if(!isfinite(ixOut) || !isfinite(iyOut)) //Probably a numerical issue
  //    return false;

  return true; //All OK
}


// Intersection line with correction of coordinates

bool PointsLineIntersectFinal( GridNode &tubeC, float x1, float x2, float y1, float y2) //Output 
{
    //http://mathworld.wolfram.com/Line-LineIntersection.html

  TVector3 dir = tubeC.m_WireDirection;

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

      tubeC.m_xDet =  tubeC.m_x;// xnom / denom;   
      tubeC.m_yDet = tubeC.m_y;//;ynom / denom;
      tubeC.m_zDet = 0.0; //z3 + 2*tubeC.m_halfLength*dir[2]*(tubeC.m_yDet  - y3)/(2*tubeC.m_halfLength * dir[1]);
      return false;
    }  else {
    distx3 = sqrt(pow(xnom/denom-x3,2) + pow(ynom/denom -y3,2));
    distx4 = sqrt(pow(xnom/denom-x4,2) + pow(ynom/denom -y4,2));
    disx3x4 = sqrt(pow(x4-x3,2) + pow(y4 -y3,2));
    if(distx3 +distx4 > disx3x4+1){
      if(distx3<distx4){
	tubeC.m_xDet = x3;   
	tubeC.m_yDet = y3;
	tubeC.m_zDet = z3 + 2*tubeC.m_halfLength*dir[2]*(tubeC.m_yDet  - y3)/(2*tubeC.m_halfLength * dir[1]);
      } else if (distx4<distx3){
	tubeC.m_xDet = x4;   
	tubeC.m_yDet = y4;
	tubeC.m_zDet = z3 + 2*tubeC.m_halfLength*dir[2]*(tubeC.m_yDet  - y3)/(2*tubeC.m_halfLength * dir[1]);
      }
      // error("Point is not on the line", xnom / denom, x3, x4, ynom / denom, y3, y4);
      return true;
    }
  }
  tubeC.m_xDet = xnom / denom;   
  tubeC.m_yDet = ynom / denom;
  tubeC.m_zDet = z3 + 2*tubeC.m_halfLength*dir[2]*(tubeC.m_yDet  - y3)/(2*tubeC.m_halfLength * dir[1]);
  //if(!isfinite(ixOut) || !isfinite(iyOut)) //Probably a numerical issue
  //    return false;

  return true; //All OK
}



// Correction intersection point 
void TubeIntersectionPointCoord(GridNode &tubeA, GridNode &tubeB)
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

  tubeA.m_xDet = ((1.0-sc) * startTubeA_x + sc * endTubeA_x);
  tubeA.m_yDet  = (1.0-sc) * startTubeA_y + sc * endTubeA_y;
  tubeA.m_zDet = (1.0-sc) * startTubeA_z + sc * endTubeA_z;

  tubeB.m_xDet = ((1.0-tc) * startTubeB_x + tc * endTubeB_x);
  tubeB.m_yDet  = (1.0-tc) * startTubeB_y + tc * endTubeB_y;
  tubeB.m_zDet = (1.0-tc) * startTubeB_z + tc * endTubeB_z;

  ////////////
  return;
}





void fit_circle(std::vector<point3D> const &pnts, CurvatureParameters &curvature)
{

  double *datax = (double*) malloc(pnts.size()*sizeof(double));
  double *datay = (double*) malloc(pnts.size()*sizeof(double));

  double xv = 0, yv = 0, r = 5;
  for(int i = 0; i < (int) pnts.size(); i++){
    datax[i] = pnts[i].m_x;
    datay[i] = pnts[i].m_y;
    xv += datax[i];
    yv += datay[i];
  }
  
  int size = (int) pnts.size();
  xv /= (double) size;
  yv /= (double) size;
  CircleData circleD(size, datax, datay);
  Circle cirini (xv,yv,r);
  Circle circle =CircleFitByHyper (circleD);
  // printf("%f, %f\n", circle.a, circle.b);
  // Circle circle;
  // CircleFitByLevenbergMarquardtFull (circleD, cirini, 0.001, circle);
  curvature.m_a = circle.a;
  curvature.m_b = circle.b;
  curvature.m_r = circle.r;
  curvature.m_E = circle.s;

  
  free(datax);
  free(datay);
}

