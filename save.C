for(unsigned int n = 0; n < nactive_queue; ++n) {
      std::vector<int> samelayer;
      std::vector<int> nextlayer;
      std::vector<int> prevlayer;

      std::vector<int> next_virt;

      std::vector<int> prev;
      curid =activeId[n];
      orid = curid;
      curindex = gr.Find(curid);
      current_Node = &Ingrid[curindex];

      n_neighbors = current_Node->m_neighbors.size();
      if(n_neighbors == 1 && current_Node->visited == 0){
	printf("\n\n\n Cur node %d, neighbors %d, visiited %d \n", curid, n_neighbors,current_Node->visited );

	PathCandidate *cand = new PathCandidate();// Create a new candidate
	cand->m_id = candidateId++;// Set id
	(cand->m_memberIdSet)->insert(curid);
	cand->m_isValid = true;
	current_Node->m_cm.push_back(cand->m_id);
	current_Node->visited = 5;

	//	neigh_ID = current_Node->m_neighbors[0];
	//next.push_back(neigh_ID);
	next.insert(next.end(),  current_Node->m_neighbors.begin(),  current_Node->m_neighbors.end());

	prev.push_back(curid);
	n_neighbors = next.size();
	bool cond = true;
	while(true){
	  if(n_neighbors == 1){

	    neigh_ID = next[0];
	    printf("1 Neighbor: node %d \n", neigh_ID);

	    neigh_index = gr.Find(neigh_ID);
	    neighnode = &Ingrid[neigh_index];
	    (cand->m_memberIdSet)->insert(neighnode->m_detID);
	    neighnode->m_cm.push_back(cand->m_id);
	    neighnode->visited = 5;

	    for(int i = 0; i < prev.size(); i ++){
	      (neighnode->m_neighbors).erase(std::remove((neighnode->m_neighbors).begin(), (neighnode->m_neighbors).end(), prev[i]), (neighnode->m_neighbors).end());
	      info("%d removed from neighbor list of %d", prev[i], neigh_ID);
	      if(i == 0)
		neighnode->parent = curid;
	    }


	    curid = neigh_ID;
	    curindex = neigh_index;
	    current_Node = neighnode;
	    next.clear();
	    prev.clear();
	    next.insert(next.end(),  current_Node->m_neighbors.begin(),  current_Node->m_neighbors.end());
	    prev.push_back(curid);
	    n_neighbors = next.size();

	  } else{

	    neigh_ID1 = next[0];
	    neigh_index1 = gr.Find(neigh_ID1);
	    neighnode1 = &Ingrid[neigh_index1];
	    
	    neigh_ID2 = next[1];
	    neigh_index2 = gr.Find(neigh_ID2);
	    neighnode2 = &Ingrid[neigh_index2];
	    printf("2 neighbors, Are these neighbors? %d, %d \n", neigh_ID1, neigh_ID2);

	    //	n_neighbors = current_Node->m_neighbors.size();
	//	for  ( int i = 0; i < current_Node.m_neighbors.size(); i++)
	//	  printf("NEIGH: %d \n",  current_Node.m_neighbors[i]);
	    //	bool cond = n_neighbors <= 2? true: false;

	      
	    if(neighnode1->IsNeighboring(neigh_ID2) || (neighnode1->m_type == GridNode::VIRTUAL_NODE && neighnode2->m_type == GridNode::VIRTUAL_NODE)){

	      if((neighnode1->m_type == GridNode::VIRTUAL_NODE && neighnode2->m_type == GridNode::VIRTUAL_NODE)){
		next_virt.clear();
		std::set_union(neighnode1->m_neighbors.begin(), neighnode1->m_neighbors.end(),
			       neighnode2->m_neighbors.begin(), neighnode2->m_neighbors.end(),
			       std::back_inserter(next_virt));
		for(int i = 0; i< prev.size(); i++)
		  (next_virt).erase(std::remove((next_virt).begin(), (next_virt).end(), prev[i]), (next_virt).end());

		if(next_virt.size() > 2)
		  break;
		else if (next_virt.size() == 2){
		  neigh_ID = next_virt[0];
		  neigh_index = gr.Find(neigh_ID);
		  neighnode =  &Ingrid[neigh_index];
		  if(!(neighnode->IsNeighboring(next_virt[1])))
		    break;
		}
	      }
	      printf("Yes, looking for common neighbors %d, %d \n", neigh_ID1, neigh_ID2);
	      (cand->m_memberIdSet)->insert(neighnode1->m_detID);
	      neighnode1->visited = 5;
	      neighnode1->m_cm.push_back(cand->m_id);
	      
	      (cand->m_memberIdSet)->insert(neighnode2->m_detID);
	      neighnode2->visited = 5;
	      neighnode2->m_cm.push_back(cand->m_id);

	      for(int i = 0; i< prev.size(); i++){
		(neighnode1->m_neighbors).erase(std::remove((neighnode1->m_neighbors).begin(), (neighnode1->m_neighbors).end(), prev[i]), (neighnode1->m_neighbors).end());
		(neighnode2->m_neighbors).erase(std::remove((neighnode2->m_neighbors).begin(), (neighnode2->m_neighbors).end(), prev[i]), (neighnode2->m_neighbors).end());
		info("%d removed from neighbor list of %d and %d", prev[i], neigh_ID1, neigh_ID2);

		
	      }
	      (neighnode1->m_neighbors).erase(std::remove((neighnode1->m_neighbors).begin(), (neighnode1->m_neighbors).end(), neigh_ID2), (neighnode1->m_neighbors).end());
	      (neighnode2->m_neighbors).erase(std::remove((neighnode2->m_neighbors).begin(), (neighnode2->m_neighbors).end(), neigh_ID1), (neighnode2->m_neighbors).end());
	      neighnode1->parent = curid;
	      neighnode2->parent = curid;

	      
	      next.clear();
	      prev.clear();
	      std::set_union(neighnode1->m_neighbors.begin(), neighnode1->m_neighbors.end(),
	      		     neighnode2->m_neighbors.begin(), neighnode2->m_neighbors.end(),
			     std::back_inserter(next));  
	      // next.insert(next.end(),  neighnode1->m_neighbors.begin(),  neighnode1->m_neighbors.end());
	      // next.insert(next.end(),  neighnode2->m_neighbors.begin(),  neighnode2->m_neighbors.end());
	      prev.push_back(neigh_ID1);
	      prev.push_back(neigh_ID2);



	      n_neighbors = next.size();
	      /*	      if(next.size() == 0){
		//	info("All neighbors, continuing");
		info("Check number of neighbors %d", next.size());
		n_neighbors = next.size();*/
	    } 

	      

	    else 
	      break;

	  }
	    for (int i = 0; i < next.size(); i++){
	      info("next on the list %d", next[i]);
	      int curidd =next[i];
	      int curindexz = gr.Find(curidd);
	      GridNode *current_Nodee = &Ingrid[curindexz];
	      for (int j = 0; j < current_Nodee->m_neighbors.size(); j++)
		info(" Has neighbors %d", current_Nodee->m_neighbors[j]);
	    }
	  // std::string dummy;
	  //  std::cout << "Enter to continue..." << std::endl;
	  //  std::getline(std::cin, dummy);
	  
	//std::vector< GridNode > possible;
	}
	  //	if (n_neighbors > 2)
	  // remainingActiveId.push_back(curid);

	//	printf(" %d, %f, %f, %f, %f, %f, %f \n", first_Node->m_detID, first_Node->m_x, first_Node->m_y, first_Node->m_z,	       first_Node->m_xDet, first_Node->m_yDet, first_Node->m_z_Det);

	temCandid.push_back(cand);

      }
      else
	remainingActiveId.push_back(curid);
	
	
      
    }








/*
	    neigh_ID1 = next[0];
	    neigh_index1 = gr.Find(neigh_ID1);
	    neighNode1 = &Ingrid[neigh_index1];
	    
	    neigh_ID2 = next[1];
	    neigh_index2 = gr.Find(neigh_ID2);
	    neighnode2 = &Ingrid[neigh_index2];
	    printf("2 neighbors, Are these neighbors? %d, %d \n", neigh_ID1, neigh_ID2);

	    //	n_neighbors = currentNode->m_neighbors.size();
	//	for  ( int i = 0; i < currentNode.m_neighbors.size(); i++)
	//	  printf("NEIGH: %d \n",  currentNode.m_neighbors[i]);
	    //	bool cond = n_neighbors <= 2? true: false;

	      
	    if(neighnode1->IsNeighboring(neigh_ID2) || (neighnode1->m_type == GridNode::VIRTUAL_NODE && neighnode2->m_type == GridNode::VIRTUAL_NODE)){

	      if((neighnode1->m_type == GridNode::VIRTUAL_NODE && neighnode2->m_type == GridNode::VIRTUAL_NODE)){
		next_virt.clear();
		std::set_union(neighnode1->m_neighbors.begin(), neighnode1->m_neighbors.end(),
			       neighnode2->m_neighbors.begin(), neighnode2->m_neighbors.end(),
			       std::back_inserter(next_virt));
		for(int i = 0; i< prev.size(); i++)
		  (next_virt).erase(std::remove((next_virt).begin(), (next_virt).end(), prev[i]), (next_virt).end());

		if(next_virt.size() > 2)
		  break;
		else if (next_virt.size() == 2){
		  neigh_ID = next_virt[0];
		  neigh_index = gr.Find(neigh_ID);
		  neighnode =  &Ingrid[neigh_index];
		  if(!(neighnode->IsNeighboring(next_virt[1])))
		    break;
		}
	      }
	      printf("Yes, looking for common neighbors %d, %d \n", neigh_ID1, neigh_ID2);
	      (cand->m_memberIdSet)->insert(neighnode1->m_detID);
	      neighnode1->visited = 5;
	      neighnode1->m_cm.push_back(cand->m_id);
	      
	      (cand->m_memberIdSet)->insert(neighnode2->m_detID);
	      neighnode2->visited = 5;
	      neighnode2->m_cm.push_back(cand->m_id);

	      for(int i = 0; i< prev.size(); i++){
		(neighnode1->m_neighbors).erase(std::remove((neighnode1->m_neighbors).begin(), (neighnode1->m_neighbors).end(), prev[i]), (neighnode1->m_neighbors).end());
		(neighnode2->m_neighbors).erase(std::remove((neighnode2->m_neighbors).begin(), (neighnode2->m_neighbors).end(), prev[i]), (neighnode2->m_neighbors).end());
		info("%d removed from neighbor list of %d and %d", prev[i], neigh_ID1, neigh_ID2);

		
	      }
	      (neighnode1->m_neighbors).erase(std::remove((neighnode1->m_neighbors).begin(), (neighnode1->m_neighbors).end(), neigh_ID2), (neighnode1->m_neighbors).end());
	      (neighnode2->m_neighbors).erase(std::remove((neighnode2->m_neighbors).begin(), (neighnode2->m_neighbors).end(), neigh_ID1), (neighnode2->m_neighbors).end());
	      neighnode1->parent = curId;
	      neighnode2->parent = curId;

	      
	      next.clear();
	      prev.clear();
	      std::set_union(neighnode1->m_neighbors.begin(), neighnode1->m_neighbors.end(),
	      		     neighnode2->m_neighbors.begin(), neighnode2->m_neighbors.end(),
			     std::back_inserter(next));  
	      // next.insert(next.end(),  neighnode1->m_neighbors.begin(),  neighnode1->m_neighbors.end());
	      // next.insert(next.end(),  neighnode2->m_neighbors.begin(),  neighnode2->m_neighbors.end());
	      prev.push_back(neigh_ID1);
	      prev.push_back(neigh_ID2);



	      n_neighbors = next.size();

	    } 
	    */
	      

	    //	    else 
	    //     break;

	  
	  /*	    for (int i = 0; i < next.size(); i++){
	      info("next on the list %d", next[i]);
	      int curIdd =next[i];
	      int curIdxz = gr.Find(curIdd);
	      GridNode *currentNodee = &Ingrid[curIdxz];
	      for (int j = 0; j < currentNodee->m_neighbors.size(); j++)
		info(" Has neighbors %d", currentNodee->m_neighbors[j]);
		}*/
	  // std::string dummy;
	  //  std::cout << "Enter to continue..." << std::endl;
	  //  std::getline(std::cin, dummy);
	  
	//std::vector< GridNode > possible;









	  /* info("Looking for a better starting point ?");
	  curId = neighId1;
	  curIdx = neighIdx1;
	  currentNode = neighNode1;
	  bool searching = true;
	  while(searching){
	    for(int i = 0; i < currentNode->m_neighbors.size(); i++){
	      int neighId = currentNode->m_neighbors[i];
	      if(neighId == prevId) continue;
	      int neighIdx  = gr.Find(neighId);
	      GridNode *neighNode = &Ingrid[neighIdx];
		  
	      if(neighNode->m_type == GridNode::VIRTUAL_NODE){
		neighId    = neighNode->m_neighbors[0] == curId ? neighNode->m_neighbors[1]: neighNode->m_neighbors[0];
		neighIdx   = gr.Find(neighId);
		neighNode  = &Ingrid[neighIdx];
		debug("This one was virtual, finding the one behind");
	      }	 
		
	      nextIds.push_back(neighId);
	    }
	    
	    if(nextIds.size() == 1){
	      info("One neighbor %d, we can continue", nextIds[0]);
	      curId = nextIds[0];
	      curIdx = gr.Find(curId);
	      currentNode = &Ingrid[curIdx];
	    } else if(nextIds.size() == 2) {
	      int testId1 = nextIds[0];
	      int testIdx1 = gr.Find(testId1);
	      GridNode* testNode1 = &Ingrid[textIdx1];

	      int testId2 = nextIds[1];
	      int testIdx2 = gr.Find(testId2);
	      GridNode* testNode2 = &Ingrid[textIdx2];
	      
	      }*/



       





	  /********** LALLALALALALALALAL */
	  
	 
	  
	  /*	  for(int i = 0; i < currentNode.size(); i++){
		  int neighId = currentNode->m_neighbors[i];
		  if(prevId != neighId){
	      

		
		  for (int i = 0; i < 2; i++){
		  int neigh_ID = first_Node->m_neighbors[i];
		  int neigh_index = gr.Find(neigh_ID);

		  GridNode *neighnode = &Ingrid[neigh_index];
		  printf("Neigh node %d \n", neigh_ID);

		  (cand->m_memberIdSet)->insert(neighnode->m_detID);
		  neighnode->m_cm.push_back(cand->m_id);
		  neighnode->visited = 5;
		  (neighnode->m_neighbors).erase(std::remove((neighnode->m_neighbors).begin(), (neighnode->m_neighbors).end(), orid), (neighnode->m_neighbors).end());


		  neighnode->parent = orid;
		  if (i == 0)
		  currentNode->visited = 5;

		  curId = neigh_ID;
		  curIdx = neigh_index;
		  currentNode = neighnode;
		  n_neighbors = currentNode->m_neighbors.size();
		  for  ( int i = 0; i < currentNode->m_neighbors.size(); i++)
	    	  printf("NEIGH: %d \n",  currentNode->m_neighbors[i]);
		  while (n_neighbors == 1){
		  neigh_ID = currentNode->m_neighbors[0];
		  neigh_index = gr.Find(neigh_ID);
		  neighnode = &Ingrid[neigh_index];
		  (cand->m_memberIdSet)->insert(neighnode->m_detID);
		  neighnode->m_cm.push_back(cand->m_id);
		  neighnode->visited = 5;

		  (neighnode->m_neighbors).erase(std::remove((neighnode->m_neighbors).begin(), (neighnode->m_neighbors).end(), curId), (neighnode->m_neighbors).end());
		  neighnode->parent = curId;

		  // currentNode->visited = 5;
		  curId = neigh_ID;
		  curIdx = neigh_index;
		  currentNode = neighnode;
		  n_neighbors = currentNode->m_neighbors.size();
		  printf("while: Neigh node %d, neighbors %d \n", curId, n_neighbors);
		  //std::vector< GridNode > possible;
	      //  std::string dummy;
	      //  std::cout << "Enter to continue..." << std::endl;
	      // std::getline(std::cin, dummy);
	
	      }*/







    /*  for(unsigned int n = 0; n < remainingActiveId.size(); ++n) {
      int curId =remainingActiveId[n];
      int curIdx = gr.Find(curId);
      GridNode *currentNode = &Ingrid[curIdx];

      std::vector< int > toVisit;
      std::vector< int > possiCand;
      std::vector< int > cluster;

      if(!currentNode->visited){
	info("Remaining node %d", curId);
	toVisit.push_back(curId);
	cluster.push_back(curId);
	currentNode->visited = 1;	      

	while(toVisit.size() != 0){
	  curId = toVisit[0];
	  curIdx = gr.Find(curId);

	  info("Current Node %d", curId);
	  toVisit.erase(toVisit.begin());
	  currentNode = &Ingrid[curIdx];
	  int n_neighbors = currentNode->m_neighbors.size();
	  for (int i = 0; i < n_neighbors; i++){
	    int neigh_ID = currentNode->m_neighbors[i];
	    // info("Neighbor %d", neigh_ID);

	    int neigh_index = gr.Find(neigh_ID);
	    GridNode *neighnode = &Ingrid[neigh_index];
	    if (neighnode->visited == 5){
	      if (!(std::find(possiCand.begin(), possiCand.end(),neighnode->m_cm[0])!=possiCand.end())){
		possiCand.push_back(neighnode->m_cm[0]);
		info("Neigh %d, One possible CM is %d",neigh_ID, neighnode->m_cm[0]);
	      }
	    }
	    else if (neighnode->visited == 2){
	      info("Neighbor %d visited", neigh_ID);
	      continue;
	    }
	    else{
	      if(neighnode->visited == 0){
		toVisit.push_back(neigh_ID);
		cluster.push_back(neigh_ID);
		info("added node %d", neigh_ID);
		neighnode->visited = 1;
	      }
	    }
	    //std::string dummy;
	    //std::cout << "Enter to continue..." << std::endl;
	    //std::getline(std::cin, dummy);
	  }
	  currentNode->visited = 2;	      
	}
	info("These nodes are neighbors to %d trackets", possiCand.size());
	for (int i = 0; i < cluster.size(); i++)
	  info("%d", cluster[i]);
	cluster.clear();
	possiCand.clear();
      }
	      
	  
      }*/

      


/* FIND REMAINING TRACKLETS */


    
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
		  cand->m_tailNode = curId;
		  for (int i = 0; i < sameLayer.size(); i++)
		    cand->m_tailNeigh.push_back(sameLayer[i]);
	      
		  for (int i = 0; i < nextLayer.size(); i++)
		    cand->m_tailNeigh.push_back(nextLayer[i]);

		  for (int i = 0; i < prevLayer.size(); i++)
		    cand->m_tailNeigh.push_back(prevLayer[i]);
		}
		else {
		  for (int i = 0; i < sameLayer.size(); i++)
		    cand->m_headNeigh.push_back(sameLayer[i]);
	      
		  for (int i = 0; i < nextLayer.size(); i++)
		    cand->m_headNeigh.push_back(nextLayer[i]);

		  for (int i = 0; i < prevLayer.size(); i++)
		    cand->m_headNeigh.push_back(prevLayer[i]);
		  cand->m_headNode = curId;
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
		      cand->m_headNeigh.push_back(candId);
		    else
		      cand->m_headNeigh.push_back(candId);

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
		    cand->m_finished = 3;
		  } else if(abs(currentNode->m_SectorLimit) > 0 || cand->m_isOnSectorLimit){
		 
		    cand->m_isOnSectorLimit= true;
		    info("Track is on sector limit, might have a connection somewhere else");
		    cand->m_finished = 1;

		  } else {
		    int firstId = cand->m_tailNode;
		    int firstIdx = gr.Find(firstId);
		    GridNode *firstNode = &Ingrid[firstIdx];
		    if(firstNode->m_LayerLimit == 1 && currentNode->m_LayerLimit == 1) {
		      info("starting and ending on same layer");
		      cand->m_finished = 3;

		    } else {
		      info("Candidate has no more neighbors, but doesn't seem finished");
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
		  cand->m_tailNode = curId;
		  for (int i = 0; i < sameLayer.size(); i++)
		    cand->m_tailNeigh.push_back(sameLayer[i]);
	      
		  for (int i = 0; i < nextLayer.size(); i++)
		    cand->m_tailNeigh.push_back(nextLayer[i]);

		  for (int i = 0; i < prevLayer.size(); i++)
		    cand->m_tailNeigh.push_back(prevLayer[i]);
		}
		else {
		  for (int i = 0; i < sameLayer.size(); i++)
		    cand->m_headNeigh.push_back(sameLayer[i]);
	      
		  for (int i = 0; i < nextLayer.size(); i++)
		    cand->m_headNeigh.push_back(nextLayer[i]);

		  for (int i = 0; i < prevLayer.size(); i++)
		    cand->m_headNeigh.push_back(prevLayer[i]);
		  cand->m_headNode = curId;
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


