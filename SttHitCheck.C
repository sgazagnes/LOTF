/*************************************
 * Author: M. Babai (M.Babai@rug.nl) *
 * Version:                          *
 * License:                          *
 *************************************/
#include <iostream>
#include <vector>

#include "CollectSttMvdPoints.h"
#include "hitcoordinate.h"
#include "gridNode.h"

#include "TFile.h"
#include "TH1.h"

void SttHitCheck(std::string const &outFileName, std::string const &dummyOutput)
{
  std::cout << "<INFO> Creating Stt tube hit histogram.\n";
  std::vector < GridNode > detNodes;
  TFile dummy(dummyOutput.c_str(), "RECREATE","Dummy outFile Created by SttHitCheck.", 9);

  std::vector < std::vector<HitCoordinate*>* >* Hit_coords = 
    CollectSttMvdPoints(detNodes, dummy);
  
  dummy.Close();
  
  TFile outFile(outFileName.c_str(),"RECREATE","Outputfile Created by SttHitCheck.", 9);

  TH1I hitHist("sttHitHist", "SttHit histogram", 4545, 0, 4545);

  for(size_t i = 0; i < Hit_coords->size(); ++i) {
    std::vector<HitCoordinate*>* currentEvtData = Hit_coords->at(i);
    for(size_t j = 0; j < currentEvtData->size(); ++j) {
      HitCoordinate* hit = currentEvtData->at(j);
      if(hit->type == HitCoordinate::STT_TYPE) {
	int DetId = hit->m_detID;
	hitHist.Fill(DetId);
      }
    }
  }
  hitHist.Write();
  // Cleanup memory
  for(size_t i = 0; i < Hit_coords->size(); ++i) {
    std::vector<HitCoordinate*>* currentEvtData = Hit_coords->at(i);
    for(size_t j = 0; j < currentEvtData->size(); ++j) {
      delete currentEvtData->at(j);
    }
    delete Hit_coords->at(i);
  }
  delete Hit_coords;

  outFile.Close();
}
