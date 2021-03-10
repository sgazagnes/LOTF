/*************************************
 * Author: M. Babai (M.Babai@rug.nl) *
 * Version:                          *
 * License:                          *
 *************************************/
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
///////////////////////
#include "TNtuple.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TH3.h"
// > 0 means create,
#define CREATE_SEPARATE_PLOTS 1

int selectCompColour(size_t compNum, size_t numColours) {
  // Set colours per component.
  int setColor = compNum % numColours;
  int color = 1;
  switch (setColor) {
  case 0:
    color = 632;
    break;
  case 1:
    color =  8;
    break;
  case 2:
    color = 600;
    break;
  case 3:
    color = 805;
    break;
  case 4:
    color = 616;
    break;
  case 5:
    color = 432;
    break;
  case 6:
    color = 400;
    break;
  case 7:
    color = 880;
    break;
    //case 6:
  default:
    color = 860;
  }
  return color;
}
/*
 * Infile contains the ntuples
 * evtNum is the event number from which we want to draw the results.
 */
void CreatePlotAllEventComponents(std::string const &InoutFile = "Tracks_output.root",
				  std::string const &dirName   = "PlotsCanvas",
				  Double_t const ww = 1100.00,
				  Double_t const hh = 1100.00,
				  size_t  const dim = 2,
				  size_t evtNum = 0)
{
  TFile inp(InoutFile.c_str(),"UPDATE");

  TNtuple* Poss   = (TNtuple*) inp.Get("Pos");
  TNtuple* MCposs = (TNtuple*) inp.Get("MCpos");
  TNtuple *Virt = (TNtuple*) inp.Get("VirtualNodesLayer");
  Virt->SetMarkerSize(0.3);

  Poss->AddFriend("MCpos");

  PosPolar->AddFriend("MCposPolar");

  // Get number of components for the current event
  
  TNtuple* grid = (TNtuple*) inp.Get("OrigGridCoord");
  int numColours = 8;
  int color = 1;
  grid->SetMarkerSize(0.3);
  // Change to directory with input events


  // Connected components
  TCanvas *c2 = new TCanvas("c2", "x-y grid with virtuals", ww, hh );
  grid->Draw("y:x","","");
  TH2F *htemp = (TH2F*)gPad->GetPrimitive("htemp");
  htemp->GetXaxis()->SetTitle("x [cm]");
  htemp->GetYaxis()->SetTitle("y [cm]");
  htemp->SetTitle("");
  c2->Update();
  // c2->SaveAs("myalgo_cm_xy.png");

 
  TCanvas *c3 = new TCanvas("c3", "z-y grid with virtuals", 300, 300 );
  grid->Draw("y:z","","");
  TH2F *htemp = (TH2F*)gPad->GetPrimitive("htemp");
  htemp->GetXaxis()->SetTitle("z [cm]");
  htemp->GetYaxis()->SetTitle("y [cm]");
  htemp->SetTitle("");
  std::stringstream EeventNumString;
  EeventNumString << evtNum;
  std::string ePlotCondition;
  for(size_t j = 0; j < nComp; j++) {
    std::stringstream ecomponentNumber;
    ecomponentNumber << j;
    // Set marker and colours per component.
    marker = 20;//(j % 34);
    int idmctr = idMatch[j];
    //printf("Id match MC %d \n", idmctr);
    if(idmctr != -1)
      color = selectCompColour(idmctr, numColours);
    else{
      color = 1;
      marker = badmarker++;
    }
    ConComps->SetMarkerColor(color);
    //  marker = 20;//(j % 34);
    //  if( marker < 20) { marker += 20; }
    // if( marker > 34) { marker = 34; }
    ConComps->SetMarkerStyle(marker);
    ConComps->SetMarkerSize(0.8);
    // Conditions
    ePlotCondition  = "(EvtNum == "  + EeventNumString.str() + ") && ";
    ePlotCondition += "(CompNum == " + ecomponentNumber.str() + ")";
    ConComps->Draw("y_Det:z_Det", ePlotCondition.c_str(), "same");
  }
  c3->Update();
  c3->SaveAs("myalgo_cm_z.png");


}
