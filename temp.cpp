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


int fitnextId(std::vector<double> &x, std::vector<double> &y, std::vector<double> &next, GridNode *goodNode, int method){

  if(method == 0){ // linear

  }  else if(method == 1){ // quadratic
    int xDir = returnDirection(  x[x.size()-2],  x[x.size()-1]);
    int yDir = returnDirection(  y[y.size()-2],  y[y.size()-1]);

    debug("We are going in x dir %d and y dir %d", xDir, yDir);

    double minDist = std::numeric_limits<double>::max();
    int goodId = -1;
    GridNode *goodNode;
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
      GridNode &node = Ingrid[curIdx];
      double xdet = (double) node.m_xDet;
      double ydet = (double) node.m_yDet;

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
      double disttube = distanceBetweenTube(node, lastNode);
      debug("Tube Distance to point %d is %lf",curId, disttube);


      int xDir = returnDirection(  x[x.size()-1], xdet);
      int yDir = returnDirection(  y[y.size()-1], ydet);

      debug("This node is in x dir %d and y dir %d", newxdir, newydir);

		
      if(minDist > currDist){
	if(labs(newxdir - xdir) > 1 || labs(newydir - ydir) > 1){
	  debug("Let's avoid going back if we can");
	} else{
		    
	  minDist = currDist;
		  
	  goodId = curId;
	  goodNode = &node;
	}
      }
	    
		
    } // FOR RIGHT NEIGHBORS
    
  }
  info("The good id is %d, and is at distance %lf", goodId, minDist);

  return goodId;
}


