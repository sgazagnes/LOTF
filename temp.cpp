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



      
    for(unsigned int l = 0; l < tracklets.size(); l++){ // Go for each tracklet
	
      PathCandidate &curCand = *(tracklets[l]);
      dbgfit("Track %d, status %d, length %d", curCand.m_id, curCand.m_finished,curCand.m_length);
	
      if (curCand.m_finished == 3 || curCand.m_length < 5 ) continue;

      int first =  curCand.m_tailNode;
      int firstIdx = gr.Find(first);
      GridNode &firstNode = Ingrid[firstIdx];

      int last =  curCand.m_headNode;
      int lastIdx = gr.Find(last);
      GridNode &lastNode = Ingrid[lastIdx];

      
      dbgfit("Tail node %d (num neigh %d),  head node %d (num neigh %d)", first, curCand.m_tailNeigh.size(), last, curCand.m_headNeigh.size());

      
      for(int k = 0; k < 2; k++){ // First Tail, then Head
	int prevId = k == 1? curCand.m_headNode: curCand.m_tailNode;
	//int curLayer = k == 1? lastNode.m_Layer: firstNode.m_Layer;
	// int layerCurDiff;
	GridNode *prevNode = k == 1? &lastNode: &firstNode;
	
	std::vector<int> *curNeigh = k == 1? &(curCand.m_headNeigh): &(curCand.m_tailNeigh);
	std::vector<unsigned int> *curMerge = k == 1? &(curCand.m_toMergeHead):&(curCand.m_toMergeTail);
	std::vector<int> next;

	if(curNeigh->size() == 0 && curMerge->size() > 0)
	  continue;


	if (curCand.m_finished == 2){
	  //  continue;
	  dbgfit("There is no neighbors, can we find some further");
	  dbgfit("Looking into %lu remaining nodes", activeId.size());

	  for(unsigned int n = 0; n < activeId.size(); ++n) {
	    int testID 		= activeId[n];
	    int testIdx = gr.Find(testID);
	    GridNode &testNode = Ingrid[testIdx];
	    float currDist =sqrt((prevNode->m_x - testNode.m_x) * (prevNode->m_x - testNode.m_x) +
				     (prevNode->m_y - testNode.m_y) * (prevNode->m_y - testNode.m_y));

	    if(currDist<5.){
		curNeigh->push_back(testID);
		dbgfit("Adding node %d to potential neighbors, into k %d", testID, k);	     
	    }
	  }

	  dbgfit("Looking into %lu  tracklets", tracklets.size());

	  for(unsigned int n = 0; n < tracklets.size(); ++n) {
	    PathCandidate &testCand = *(tracklets[n]);
	
	    if (testCand.m_finished == 3 || n == l) continue;

	    int tailId =  testCand.m_tailNode;
	    int  tailIdx = gr.Find(tailId);
	    GridNode &tailNode = Ingrid[tailIdx];

	    int headId =  testCand.m_headNode;
	    int headIdx = gr.Find(headId);
	    GridNode &headNode = Ingrid[headIdx];

	    GridNode Dummy;
	    
	    double currDistTail = IntersectionPointSkePar(gr, *prevNode, tailNode, Dummy);//sqrt((prevNode->m_x - tailNode.m_x) * (prevNode->m_x - tailNode.m_x) +
				  //   (prevNode->m_y - tailNode.m_y) * (prevNode->m_y - tailNode.m_y));
	    double currDistHead = IntersectionPointSkePar(gr, *prevNode, headNode, Dummy);//sqrt((prevNode->m_x - headNode.m_x) * (prevNode->m_x - headNode.m_x) +
	      //(prevNode->m_y - headNode.m_y) * (prevNode->m_y - headNode.m_y));
	    dbgfit("%d %f, %d %f", tailId, currDistTail, headId, currDistHead);
	    if(currDistTail < 5. || currDistHead < 5.){
	      if(currDistTail <= currDistHead){
		curNeigh->push_back(tailId);
		dbgfit("Adding node %d to potential neighborsl", tailId);
	      } else {
		curNeigh->push_back(headId);
		dbgfit("Adding node %d to potential neighbors of head", headId);
	      }
	    }
	  }
	  
	}

	
	if(curNeigh->size() == 0){
	  dbgfit("No good candidate has been found");
	  continue;
	}
	


	k == 1? dbgfit("HEAD : Fitting next neighbors "): dbgfit("TAIL : Fitting next neighbors");
	      
	for(size_t i = 0; i  < curNeigh->size(); i++){ 
	  int id = curNeigh->at(i);
	  int idx = gr.Find(id);
	  GridNode &node = Ingrid[idx];
	  if(node.m_type == GridNode::VIRTUAL_NODE){
	    //	  for (size_t j = 0; j < node.m_neighbors.size(); j++) {
	       // Remove second order neighbors from virtual
	      int neigh1 = node.m_neighbors[0];
	      int neigh2 = node.m_neighbors[1];
	      curNeigh->erase(std::remove(curNeigh->begin(), curNeigh->end(), neigh1), curNeigh->end());
	      curNeigh->erase(std::remove(curNeigh->begin(), curNeigh->end(), neigh2), curNeigh->end());
	      //  }
	  }
	}
	  
	next.insert(next.end(),  curNeigh->begin(),  curNeigh->end());

	bool cond = next.size() > 0? true: false;
	bool adding = true;

	std::vector<int> virt;
	std::vector<int> *trk = curCand.m_memberList;

	int potCm = -1;
	int id = k == 1? trk->at(trk->size() - 2) : 1;
	int idx = gr.Find(id);
	GridNode &node = Ingrid[idx];		
	   

	while (cond){
	    
	  GridNode *goodNode;	    
	    
	  int goodId = fitNextId(gr, curCand, next, k);
	    
	  if (goodId == -1) {
	    dbgfit("No good candidates have been found, stop");
	    dbgfit("Current cm %d: length is %d,  tail node %d  head node %d  Min layer %d, Max layer %d. ", curCand.m_id, curCand.m_length, curCand.m_tailNode, curCand.m_headNode, curCand.m_minLayer, curCand.m_maxLayer);
	    if((curCand.m_minLayer == 0 && curCand.m_maxLayer > 21) || (firstNode.m_LayerLimit == 1 && lastNode.m_LayerLimit == 1)){		 
	      dbgfit("track goes through all layers or makes a loop, likily finished");		 
	      curCand.m_finished = 3;		 
	    } else {
	      curCand.m_finished = 2;
	      cond = false;
	      break;
	    }
	  }


	      
	  int goodIdx = gr.Find(goodId);
	  goodNode = &Ingrid[goodIdx];

	      
	 
	    
	  // Automatically add next neighbord if we just found a virtual node
	    
	  if(goodNode->m_type == GridNode::VIRTUAL_NODE){
	    curCand.insertNewNode(gr,goodNode, k == 0? curCand.m_memberList->begin(): curCand.m_memberList->end());
	    visited[goodId] = 4;
	    int neighId = curCand.isInCandidate(goodNode->m_neighbors[0]) ? goodNode->m_neighbors[1]:goodNode->m_neighbors[0];
	    int neighIdx = gr.Find(neighId);
	    goodNode = &Ingrid[neighIdx];
	    dbgfit("We added a virtual node, also adding the next neighbord %d", neighId);
	    goodId = neighId;
	  }
	    
	    
	  if(visited[goodId] == 4 && !goodNode->m_LayerLimit ){

	    potCm = goodNode->m_cm[0];
	    const auto p = std::find_if(tracklets.begin(), tracklets.end(),
					[potCm](const PathCandidate *obj){ return obj->m_id == potCm; } );

	    PathCandidate &neighCand = *(*p);

	    bool willMerge = neighCand.m_toMergeHead.size() == 0 && neighCand.m_toMergeTail.size() == 0?false:true;
	    /*     for(int i = 0; i<tracklets.size(); i++){
		   if(sayYes[potCm][i] != 0)
		   willMerge = true;
		   }*/

	    //	    if(!willMerge){

	      std::vector<int>::iterator it = std::find((neighCand.m_memberList)->begin(), (neighCand.m_memberList)->end(), goodId);
	      int index = std::distance((neighCand.m_memberList)->begin(), it);
	      int id = neighCand.m_memberList->at(index);
	      int nextNeigh;
	      int dir = 0;

	      dbgfit("This node %d already belongs to a CM (%d), tail Node %d, and head Node %d", id, goodNode->m_cm[0],neighCand.m_tailNode,neighCand.m_headNode);
	      if (std::find(curCand.m_toMergeHead.begin(), curCand.m_toMergeHead.end(), goodNode->m_cm[0]) != curCand.m_toMergeHead.end() || std::find(curCand.m_toMergeTail.begin(), curCand.m_toMergeTail.end(), goodNode->m_cm[0]) != curCand.m_toMergeTail.end() ){
		dbgfit("We were already planning to merge, then let's stop here");
		break;
	      }  

	      if(id != neighCand.m_headNode && id != neighCand.m_tailNode){
		dbgfit("The index is neither the tail or the head, we should continue");
	      }else{
		if(id == neighCand.m_headNode){
		  dbgfit("We will merge in the tail direction");
		  nextNeigh = index-1;
		  dir = k == 1? 3: 1; // head to head or tail to head
		  k == 1? curCand.m_toMergeHead.push_back(potCm):  curCand.m_toMergeTail.push_back(potCm);
		  /*  if((neighCand.m_toMergeHead).size() > 0){
		    dbgfit("WE WERE ALREADY MERGING IN THAT DIRECTION ?");
		    break;
		    }*/
		  neighCand.m_toMergeHead.push_back(curCand.m_id);
		  sayYes[curCand.m_id][potCm] = dir;
		} else if ( id == neighCand.m_tailNode ) {
		  dbgfit("We will merge in the head direction");
		  nextNeigh= index+1;
		  dir = k == 1? 2: 0; // head to tail or tail to tail
		  k == 1? curCand.m_toMergeHead.push_back(potCm):  curCand.m_toMergeTail.push_back(potCm);
		  /*  if((neighCand.m_toMergeTail).size() > 0){
		    dbgfit("WE WERE ALREADY MERGING IN THAT DIRECTION ?");
		    break;
		    }*/
		  neighCand.m_toMergeTail.push_back(curCand.m_id);
		  sayYes[curCand.m_id][potCm] = dir;
		}
		  
		float angle_r = returnAngle(prevNode->m_r, goodNode->m_r, neighCand.m_r[nextNeigh], (prevNode->m_thetaDeg+180.)/360., (goodNode->m_thetaDeg+180.)/360., (neighCand.m_theta[nextNeigh]+180.)/360.);
		  
		float angle_xy = returnAngle(prevNode->m_x, goodNode->m_x, neighCand.m_x[nextNeigh], prevNode->m_y, goodNode->m_y, neighCand.m_y[nextNeigh]);
		  
		dbgfit("Angle r track %f", angle_r);
		dbgfit("Angle xy track %f", angle_xy);
		  
		break;
	      }
	      //    // } else 
	      //     dbgfit("we did not merge because we were in the middle of an other track");
	  } // Visited good id 4
	  /*else {
	    debug("we must find which is the direction");
	    int currr =  neighCand.m_r[index-1] - goodNode->m_r;
	    int prevr = goodNode->m_r - prevNode->m_r;
	    if((currr >= 0 && prevr >= 0) || (currr < 0 && prevr < 0)){
	    debug("Middle node, We are checking in the tail direction");
	    nextNeigh= index-1;
	    dir = k == 1? 3: 1; // head to head or tail to head
	    } else {
	    debug("Middle node, We are checking in the head direction");
	    nextNeigh= index+1;
	    dir = k == 1? 2: 0; // head to head or tail to head
	    }
	    }*/



	       
	  /*  if(potCm == goodNode->m_cm[0]){
		  
	      debug("Two nodes in a row belong to the same cm %d, let's stop and merge them later", potCm);
	      cond = false;
	      const auto p = std::find_if(tracklets.begin(), tracklets.end(),
	      [potCm](const PathCandidate *obj){ return obj->m_id == potCm; } );
		

	      curCand.m_toMergeHead.push_back(potCm);

	      if(prevId == (*p)->m_headNode){
	      (*p)->m_toMergeHead.push_back(curCand.m_id);
	      debug("Merging head to head");
	      }
		  
	      else if(prevId ==  (*p)->m_tailNode){
	      (*p)->m_toMergeTail.push_back(curCand.m_id);
	      debug("Merging head to tail");
	      }

	      else {
	      std::vector<int>::iterator it = std::find(((*p)->m_memberList)->begin(), ((*p)->m_memberList)->end(), prevId);
	      int index = std::distance(((*p)->m_memberList)->begin(), it);
	      if(index >= (tracklets[potCm]->m_memberList)->size()/2){
	      debug("Merging head to head (to check)");
	      ( (*p)->m_toMergeHead).push_back(curCand.m_id) ;
	      } else{
	      debug("Merging head to tail (to check)");
	      ( (*p)->m_toMergeTail).push_back(curCand.m_id);
	      }
	      }
	      break; 
	      //visitedTracks[curCand.m_id] = 1;
	      // visitedTracks[potCm] = 1;
	      } else {
	      potCm = goodNode->m_cm[0];
	      }	*/	
	   
	  //else
	  // potCm = -1;

	  // goodIdx = gr.Find(goodId);
	  curCand.insertNewNode(gr, goodNode, k == 0? curCand.m_memberList->begin(): curCand.m_memberList->end());
	  visited[goodId] = 4;


	  //   curLayer = goodNode->m_Layer;

	  next.clear();

	  for(size_t i = 0; i < goodNode->m_neighbors.size(); i++){
	      
	    int neighId = goodNode->m_neighbors[i];
	    //  debug("Node %d has one neig %d, visited %d", goodId, neighId, visited[neighId]);

	    if(curCand.isInCandidate(neighId)) continue;

	      

	    int neighIdx = gr.Find(neighId);
	    GridNode *neighNode = &Ingrid[neighIdx];

	    next.push_back(neighId);
		
	  }// for neighbors

	  if(goodNode->parent != -1){
	    next.push_back(goodNode->parent);
	    dbgfit("Also adding this node %d", goodNode->parent );
	  }


	  if(next.size() > 0){	      
	    prevId = goodId;
	    prevNode = goodNode;
	    virt.clear();
	    if(goodNode->m_Layer == 0){
	      dbgfit("Reached end ?");
	      curCand.m_finished = 2;
	      cond = false;
	    }			 
	  }

	  else {
	      
	    dbgfit("we have no more neighbors");
	    cond = false;
	      
	    /*if(potCm != -1){
		
	      debug("Last node added belonged to a cm %d, let's see if we should merge them", potCm);
		
	      curCand.m_toMergeHead.push_back(potCm);

	      const auto p = std::find_if(tracklets.begin(), tracklets.end(), [potCm](const PathCandidate *obj){ return obj->m_id == potCm; } );
		

	      if(goodId == (*p)->m_headNode){
	      ((*p)->m_toMergeHead).push_back(curCand.m_id);
	      debug("Merging head to head");
	      }
	      else if(goodId == (*p)->m_tailNode){
	      ((*p)->m_toMergeTail).push_back(curCand.m_id);
	      debug("Merging head to tail");
	      }else{
	      std::vector<int>::iterator it = std::find(((*p)->m_memberList)->begin(), ((*p)->m_memberList)->end(), goodId);
	      int index = std::distance(((*p)->m_memberList)->begin(), it);
	      if(index >= ((*p)->m_memberList)->size()/2){
	      debug("Merging head to head (to check)");
	      ((*p)->m_toMergeHead).push_back(curCand.m_id) ;
	      } else{
	      debug("Merging head to tail (to check)");
	      ((*p)->m_toMergeTail).push_back(curCand.m_id) ;
	      }
	      }
	      //	visitedTracks[curCand.m_id] = 1;
	      //	visitedTracks[potCm] = 1;
	      } else {
	      info("No potential neighbors ?");
	      curCand.m_finished = 2;
	      }*/
		     
	  } // next,size() == 0
			  	      
	    //   cond = false;
	} // WHILE COND
      } // for k 2
      dbgfit("Finished with current tracklet, moving to next \n");
    } // for k trackets
    // timing("Fitting phase ended. Real time %f s, CPU time %f s.", timer.RealTime(),timer.CpuTime());
