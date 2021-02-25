
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

#include "phzinterp.h"
#include "gridNode.h"
#include "logc.h"
#include "simon_functions.h"



void ZCoordinates(std::vector < PathCandidate* > &tracklets)
{
 
  for(unsigned int l = 0; l < tracklets.size(); l++){
     
    PathCandidate &curCand = *(tracklets[l]);
    if(!curCand.m_isValid)
      continue;

    dbgtrkz("Re-determining the z coordinates for the track", curCand.m_id);

    std::vector<int> const *vect = curCand.m_memberList;

  
    std::vector<double> &x =  curCand.m_x;
    std::vector<double> &y =  curCand.m_y;
    std::vector<double> &z =  curCand.m_z;
    std::vector<GridNode> &anchors =  curCand.m_anchors;

    std::vector< double >  xPts;
    std::vector< double >  yPts;
    std::vector< double >  zPts;

    int lastid = -1;
    int firstid = -1;
    for(size_t i = 0; i < anchors.size(); i++){
      GridNode  &node = anchors[i];

      if(node.m_weight && node.m_z_Det != 0.0){
	xPts.push_back((double) node.m_xDet);
	yPts.push_back((double) node.m_yDet);
	double zadd = (double) node.m_z_Det;
	if(zadd > 105)
	  zadd = 100;
	else if(zadd < -35)
	  zadd = -30;
	zPts.push_back((double) zadd);
	lastid = node.m_detID;
	firstid = firstid == -1? node.m_detID: firstid;
      }
    }
    if(xPts.size()>1){
    
      dbgtrkz("Coordinates found for fitting");

      int dir = anchors[0].m_detID < anchors[anchors.size()-1].m_detID || anchors[0].m_Layer == 0 ?
				     1 : -1;
      if(dir == -1){
	xPts.push_back(0.);
	yPts.push_back(0.);
	zPts.push_back(0.);

      } else {
	xPts.insert(xPts.begin(), 0.);
	yPts.insert(yPts.begin(), 0.);
	zPts.insert(zPts.begin(), 0.);
      }

      
      // for(size_t i = 0; i < zPts.size(); i++)
      //	dbgtrkz("x %lf, y %lf, z %lf", xPts[i],  yPts[i],  zPts[i]);


      std::vector<double> p;
      std::vector<double> distxy;

      double newval;
      double distxx = 0;
      p.push_back(0.);
      distxy.push_back(0.);
      int start =1;

      for (size_t i = 0; i <  zPts.size()-1; i++){ 
	newval = p[i] + sqrt(pow(zPts[i+1]-zPts[i],2.) + pow(yPts[i+1]-yPts[i],2.));
	p.push_back(newval);
      }
      dbgtrkz("Last id is %d, first id is %d", lastid, firstid);

     if(dir == -1){
	distxx += sqrt(x.back()*x.back()+y.back()*y.back());
	for(int i =vect->size()-1; i >= 0; i--){
	  distxx +=  sqrt(pow(x[i-1]-x[i],2.) + pow(y[i-1]-y[i],2.));
	  if(vect->at(i) == firstid)
	    break;
	}
      }else{
	distxx += sqrt(x[0]*x[0]+y[0]*y[0]);
	for(size_t i =0; i < lastid; i++){
	  if(vect->at(i) == lastid)
	    break;
	  distxx +=  sqrt(pow(x[i+1]-x[i],2.) + pow(y[i+1]-y[i],2.));
	}
      }
      double *z_coef = polyFit(p, zPts, 1);
      double *y_coef = polyFit(p, yPts, 1);

      dbgtrkz("Fitted z track is %f + p%f", z_coef[0], z_coef[1]);
      dbgtrkz("Fitted y track is %f + p%f", y_coef[0], y_coef[1]);
      //dbgtrkz("p: %lf, %lf", p[start], p[p.size()-1]);

      double pfirst = (y[0] - y_coef[0])/y_coef[1];
      double Node_distance =0;
      double z0 = dir == -1? zPts.back(): zPts[0];//z_coef[0]+z_coef[1]*p[start];
      double zf = dir == -1? z_coef[0]:z_coef[0]+z_coef[1]*p[p.size()-1];
      if(zf > 110)
	zf = 105;
      //bgtrkz("z: %lf, %lf", z0, zf);
      double y0 =y_coef[0]+y_coef[1]*p[start];
      double yf =y_coef[0]+y_coef[1]*p[p.size()-1];
      // dbgtrkz("y: %lf, %lf", y0, yf);

      double totxydist =  distxx;
      // dbgtrkz("Total fist %lf, slope %lf", totxydist, (zf-z0)/totxydist);
      bool change =0;
      if(dir == -1){
	for(int i = vect->size()-1; i >=0; i--){
	  if(i < vect->size()-1) Node_distance += sqrt(pow(x[i]-x[i+1],2)+pow(y[i]-y[i+1],2));//-xPts[start]-yPts[start]
	  else Node_distance += sqrt(pow(x[i],2)+pow(y[i],2));
	
	  z[i] = z0 + (zf-z0)*Node_distance/totxydist;
	  //	  dbgtrkz("Dist cur %lf", Node_distance);
	  //	  dbgtrkz("New %f, %f, %f", x[i], y[i], z[i]);
	}
      }else{
	for(size_t i = 0; i <vect->size(); i++){
	  if(i >= 0) Node_distance += sqrt(pow(x[i]-x[i-1],2)+pow(y[i]-y[i-1],2));//-xPts[start]-yPts[start]
	  else Node_distance += sqrt(pow(x[i],2)+pow(y[i],2));
	
	  z[i] = z0 + (zf-z0)*Node_distance/totxydist;
	  //dbgtrkz("Dist cur %lf", Node_distance);
	 // dbgtrkz("New %f, %f, %f", x[i], y[i], z[i]);
	}
      }
  
      
    } else
      dbgtrkz("No enough good coord found");

    dbgtrkz("Finished cur track \n");
  }
}
