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


void fitZCoordinates(CoordGrid &hitMap, PathCandidate *trk)
{
 
 
  dbgtrkz("Re-determining the z coordinates for the track", trk->m_id);

  // MVD list of available elements in the graph.
  // std::vector< GridNode > &MVDPoints = hitMap.m_MVD_grid;
  //
  // STT lists of elements available in the graph.
  std::vector< GridNode > &Ingrid = hitMap.m_grid;
  // Local variables
  std::set<int>::iterator stt_It;
  // Track (Connected components) )loop

  // Current track object  std::vector<double> z =  std::vector<double>( trk->m_z );

  std::vector<double> &x =  trk->m_x;
  std::vector<double> &y =  trk->m_y;
  std::vector<double> &z =  trk->m_z;
  std::vector<int> const *vect = trk->m_memberList;

  // Fnd the most outer virtual node or the max Z
  // size_t maxLayer = std::numeric_limits<size_t>::min();
  // float maxVirtZVal = std::numeric_limits<float>::min();
  //float  lastVirtualZCoord = 0.00;
  for( size_t i = 0; i < z.size(); i++){
    printf("%d, %lf \t", vect->at(i), z[i]);
  }
  printf("\n");

  std::vector<GridNode*>  virt;
  std::vector<GridNode*>  prevLayer;
  std::vector<GridNode*>  nextLayer;
  std::vector<GridNode*>  interLayer1;
  std::vector<GridNode*>  interLayer2;
  int detect = 0;
  int idxFirstVirt = 0, idxFirstSkewed;
  double xvirtfirst = 0, yvirtfirst = 0, zvirtfirst = 0;
  double xvirtsecond = 0, yvirtsecond = 0, zvirtsecond = 0;
  double lastx=0, lasty =0;
  GridNode Dummy_coord;

  for(size_t i = 0; i < vect->size(); i++){
    int nodeID = vect->at(i);
    size_t node_index = hitMap.Find(nodeID);
    GridNode  *node = &Ingrid[node_index];

    if( detect < 2 && node->m_type == GridNode::VIRTUAL_NODE) {
      virt.push_back(node);
      int curId = node->m_neighbors[0];
      size_t curIdx = hitMap.Find(curId);
      GridNode  *curNode = &Ingrid[curIdx];
      lastx = curNode->m_x;
      lasty = curNode->m_y;
      if(curNode->m_type == GridNode::STT_TYPE_SKEW){
	curId = node->m_neighbors[1];
	curIdx = hitMap.Find(curId);
	curNode = &Ingrid[curIdx];
	lastx = curNode->m_x;
	lasty = curNode->m_y;
      }
      if(detect == 0){
	idxFirstVirt = i;
	//	prevLayer.push_back(prevnode);
	//	nextLayer.push_back(nextnode);
	detect = 1;
      } 
	
    }	
    else if(detect == 1){
      
      dbgtrkz("Number of virtual nodes found %d", virt.size());


      for(size_t j = 0; j< virt.size(); j++){
	xvirtfirst += virt[j]->m_x;
	yvirtfirst += virt[j]->m_y;
	zvirtfirst += virt[j]->m_z_Det;
      }
      xvirtfirst /= virt.size();
      yvirtfirst /= virt.size();
      zvirtfirst /= virt.size();
      //xvirtfirst = (x[i] + lastx)/2.0;
      //yvirtfirst = (y[i] + lasty)/2.0;

      dbgtrkz("average coord of virt layer %lf %lf", xvirtfirst, yvirtfirst);
      for(size_t j = idxFirstVirt; j < idxFirstVirt + virt.size(); j++){
	x[j] = xvirtfirst;
	y[j] = yvirtfirst;
	z[j] = zvirtfirst;
      }
      
      virt.clear();
      // interLayer1.insert(insertLayer1.end(),  nextLayer->begin(),  nextLayer->end());
      //nextLayer.clear();
      interLayer1.push_back(node);
      idxFirstSkewed = i;
      dbgtrkz("Adding %d on interlayer 1",node->m_detID);
      detect = 2;
      
    } else if(detect == 2 &&  node->m_type != GridNode::VIRTUAL_NODE){
      if (node->m_Layer == interLayer1[0]->m_Layer && !(std::find(std::begin(interLayer1), std::end(interLayer1), node) != std::end(interLayer1))) {
	interLayer1.push_back(node);
	dbgtrkz("Adding %d on interlayer 1",node->m_detID);
      } else if (labs(node->m_Layer - interLayer1[0]->m_Layer)<2) {
	interLayer2.push_back(node);
	dbgtrkz("Adding %d on interlayer 2",node->m_detID);
      } else {
	dbgtrkz("We are missing a virtual for %d",node->m_detID);
	GridNode tubeA = *(interLayer2.back());
	GridNode tubeB = *(node);
	info("%d %d", tubeA.m_detID , tubeB.m_detID); 

	IntersectionPoint(hitMap, tubeA, tubeB, Dummy_coord);
	Dummy_coord.m_detID      = 100000;
	dbgtrkz("Dummy_coord %d",Dummy_coord.m_detID ); 
	std::pair<float, float> r_Theta;
	float theta_deg = Cartesian_To_Polar(Dummy_coord.m_xDet, Dummy_coord.m_yDet, r_Theta);
	Dummy_coord.m_r = r_Theta.first;
	Dummy_coord.m_thetaDeg = theta_deg;
	xvirtsecond = Dummy_coord.m_x;
	yvirtsecond = Dummy_coord.m_y;
	zvirtsecond = Dummy_coord.m_z;
	idxFirstVirt = -1;
	detect=3;
	i--;
	virt.push_back(&Dummy_coord);
	//	interLayer2.push_back(node);// to remove
      }
    } else if(detect == 2  && node->m_type == GridNode::VIRTUAL_NODE){
      idxFirstVirt = i;
      virt.push_back(node);
      detect = 3;
      lastx = x[i-1];
      lasty = y[i-1];
    } else if(detect == 3 && node->m_type == GridNode::VIRTUAL_NODE){
      virt.push_back(node);
    } else if(detect == 3 && node->m_type != GridNode::VIRTUAL_NODE){

      
      dbgtrkz("We have the intermediate layer. 1 %d, 2 %d, virtual further  %d", interLayer1.size(), interLayer2.size(), virt.size());

      if(idxFirstVirt != -1){
	for(size_t j = 0; j< virt.size(); j++){
	  xvirtsecond += virt[j]->m_x;
	  yvirtsecond += virt[j]->m_y;
	  zvirtsecond += virt[j]->m_z_Det;
	}
	xvirtsecond /= virt.size();
	yvirtsecond /= virt.size();
	zvirtsecond /= virt.size();
	//xvirtsecond = (x[i] + lastx)/2.0;
	//yvirtsecond = (y[i] + lasty)/2.0;

	for(size_t j = idxFirstVirt; j < idxFirstVirt + virt.size(); j++){
	  x[j] = xvirtsecond;
	  y[j] = yvirtsecond;
	  z[j] = zvirtsecond;
	}
      }
      virt.clear();


      dbgtrkz("Virt point 1 %lf %lf %lf \n virt point 2 %lf, %lf, %lf", xvirtfirst, yvirtfirst, zvirtfirst,xvirtsecond, yvirtsecond, zvirtsecond);
	
	if(interLayer1.size() <= 3 ){

	  double layerX = interLayer1[0]->m_x, layerY= interLayer1[0]->m_y, layerZ =interLayer1[0]->m_z ;
	  TVector3 dir = interLayer1[0]->m_WireDirection; 
	  double hl = interLayer1[0]->m_halfLength;
	  if(interLayer1.size() > 1){
	    for(size_t j = 1; j < interLayer1.size(); j++){
	      layerX += interLayer1[j]->m_x;
	      layerY += interLayer1[j]->m_y;
	      layerZ += interLayer1[j]->m_z;
	    }
	    layerX /= interLayer1.size();
	    layerY /= interLayer1.size();
	    layerZ /= interLayer1.size();
	  }



	  // debug(" Start x1 %lf, end x1 %lf, start y1 %lf, end y 1 %lf", layerX - hl*dir[0], layerX + hl*dir[0],  layerY - hl*dir[1],  layerY + hl*dir[1]);
	  double coef = IntersectionXY(layerX - hl*dir[0], layerX + hl*dir[0],  layerY - hl*dir[1],  layerY + hl*dir[1], xvirtfirst,xvirtsecond ,yvirtfirst, yvirtsecond);
	  double newX =  (1.0-coef) * (layerX - hl*dir[0]) + coef * (layerX + hl*dir[0]);
	  double newY =  (1.0-coef) * (layerY - hl*dir[1]) + coef * (layerY + hl*dir[1]);
	  double newZ =  (1.0-coef) * (layerZ - hl*dir[2]) + coef * (layerZ + hl*dir[2]);
      
	  dbgtrkz("Coef %lf, Intersection point on layer 1 is %lf, %lf, %lf", coef, newX, newY, newZ);

	  double diffX = newX - layerX;
	  double diffY = newY - layerY;
	  double diffZ = newZ - layerZ;
      
	  for(size_t j = idxFirstSkewed, k =0; j < idxFirstSkewed + interLayer1.size(); j++, k++){
	    x[j] = interLayer1[k]->m_x + diffX;
	    y[j] = interLayer1[k]->m_y + diffY;
	    z[j] = newZ;
	    interLayer1[k]->m_z_Det = newZ;
	    dbgtrkz("Node %d has new coord %lf, %lf %lf", vect->at(j), x[j], y[j], z[j]);
	  }
      

	  if(interLayer2.size()> 0){
      
	  layerX = interLayer2[0]->m_x;
	  layerY = interLayer2[0]->m_y;
	  layerZ = interLayer2[0]->m_z;
	  dir    = interLayer2[0]->m_WireDirection; 
	  hl     = interLayer2[0]->m_halfLength;
	  if(interLayer2.size()> 1){
	    for(size_t j = 1; j < interLayer2.size(); j++){
	      layerX += interLayer2[j]->m_x;
	      layerY += interLayer2[j]->m_y;
	      layerZ += interLayer2[j]->m_z;
	    }
	    layerX /= interLayer2.size();
	    layerY /= interLayer2.size();
	    layerZ /= interLayer2.size();
	  }
	  coef = IntersectionXY(layerX - hl*dir[0], layerX + hl*dir[0],  layerY - hl*dir[1],  layerY + hl*dir[1], xvirtfirst,xvirtsecond ,yvirtfirst, yvirtsecond);
	  newX =  (1.0-coef) * (layerX - hl*dir[0]) + coef * (layerX + hl*dir[0]);
	  newY =  (1.0-coef) * (layerY - hl*dir[1]) + coef * (layerY + hl*dir[1]);
	  newZ =  (1.0-coef) * (layerZ - hl*dir[2]) + coef * (layerZ + hl*dir[2]);
      
	  dbgtrkz("Coef %lf, Intersection point on layer 2 is %lf, %lf, %lf", coef, newX, newY, newZ);

	  diffX = newX - layerX;
	  diffY = newY - layerY;
	  diffZ = newZ - layerZ;
      
	  for(size_t j = idxFirstSkewed+interLayer1.size(), k =0; j < idxFirstSkewed +interLayer1.size()+ interLayer2.size(); j++, k++){
	    x[j] = interLayer2[k]->m_x + diffX;
	    y[j] = interLayer2[k]->m_y + diffY;
	    z[j] = newZ;
	    interLayer2[k]->m_z_Det = newZ;

	    dbgtrkz("Node %d has new coord %lf, %lf %lf", vect->at(j), x[j], y[j], z[j]);
		  
	  }
	  }
	  
      } else {
	dbgtrkz("Ambiguous, too many nodes ?");
	for(size_t j = idxFirstSkewed, k =0; j < idxFirstSkewed + interLayer1.size(); j++, k++){
	  x[j] = interLayer1[k]->m_x;
	  y[j] = interLayer1[k]->m_y;
	  z[j] = -1;
	  dbgtrkz("Node %d has new coord %lf, %lf %lf", vect->at(j), x[j], y[j], z[j]);
	}
      
      }
      
      interLayer1.clear();
      interLayer2.clear();
      xvirtfirst = xvirtsecond;
      yvirtfirst = yvirtsecond;
      xvirtsecond = yvirtsecond = 0;

      if(node->m_type != GridNode::STT_TYPE_SKEW){
	dbgtrkz("Need to interpolate now\n\n");
	break;
      }
      interLayer1.push_back(node);
      dbgtrkz("Adding %d on interlayer 1",node->m_detID);
      idxFirstSkewed = i;
      detect = 2;
    }
      
    
    // END IF virtual
  // Determine the largest Z value (Virtual node)
  }
}
