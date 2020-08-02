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

      
