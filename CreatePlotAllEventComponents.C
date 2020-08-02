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
    color = 1;
    break;
  case 1:
    color =  632;
    break;
  case 2:
    color = 418;
    break;
  case 3:
    color = 600;
    break;
  case 4:
    color = 923;
    break;
  case 5:
    color = 863;
    break;
    //case 6:
  default:
    color = 801;
  }
  return color;
}
/*
 * Infile contains the ntuples
 * evtNum is the event number from which we want to draw the results.
 */
void CreatePlotAllEventComponents(std::string const &InoutFile = "Tracks_output.root",
				  std::string const &dirName   = "PlotsCanvas",
				  Double_t const ww = 900.00,
				  Double_t const hh = 900.00,
				  size_t  const dim = 2,
				  size_t evtNum = 0)
{
  TFile inp(InoutFile.c_str(),"UPDATE");

  // Get COnnected components
  TNtuple* ConComps = (TNtuple*) inp.Get("ConnectedComponents");
  
  // Get number of components for the current event
  float numComponents;

  TTree* EvtNumComps = (TTree*) inp.Get("ComponentPerEvt");
  EvtNumComps->SetBranchAddress("numComponents", &numComponents);

  EvtNumComps->GetEntry(evtNum);
  unsigned int nComp = static_cast<unsigned int>(numComponents);

  // Fetch data and make plots
  TNtuple* grid = (TNtuple*) inp.Get("ExtendedGrid");
  int numColours = 7;
  int color = 1;
  
  // Change to directory with input events
  inp.cd("InputEvents");

  std::stringstream evN;
  std::string tupName = "Evt";
  evN << evtNum;
  tupName = tupName + "_" + evN.str() + "_CoordsTuple";
  tupName = "InputEvents/" + tupName;

  TNtuple* pos = (TNtuple*) inp.Get(tupName.c_str());

  // Create and set the canvas name and properties
  std::string canvasName = "Event_" + evN.str() + "_PlotCanvas";
  std::string descr = "All Collected Coordinates Event " + evN.str();

  TCanvas *c1 = new TCanvas(canvasName.c_str(), descr.c_str(), ww, hh );

  int marker = 0;
  
  //_____________ BEGIN TEMPORARY PLOTS May be DELETED
#if(CREATE_SEPARATE_PLOTS > 0)
  {
    // Connected components
    TCanvas *c2 = new TCanvas("c2", "Connected Components", 300, 300 );
    grid->Draw("y:x","","");
    TH2F *htemp = (TH2F*)gPad->GetPrimitive("htemp");
    htemp->GetXaxis()->SetTitle("x [cm]");
    htemp->GetYaxis()->SetTitle("y [cm]");
    htemp->SetTitle("");
    std::stringstream EeventNumString;
    EeventNumString << evtNum;
    std::string ePlotCondition;
    for(size_t j = 0; j < nComp; j++) {
      std::stringstream ecomponentNumber;
      ecomponentNumber << j;
      // Set marker and colours per component.
      color = selectCompColour(j, numColours);
      ConComps->SetMarkerColor(color);
      marker = (j % 34);
      if( marker < 20) { marker += 20; }
      if( marker > 34) { marker = 34; }
      ConComps->SetMarkerStyle(marker);
      ConComps->SetMarkerSize(0.8);
      // Conditions
      ePlotCondition  = "(EvtNum == "  + EeventNumString.str() + ") && ";
      ePlotCondition += "(CompNum == " + ecomponentNumber.str() + ")";
      ConComps->Draw("y:x", ePlotCondition.c_str(), "same");
    }
    c2->Update();
  }
  // Z-reconstructed Connected components
  {
    TCanvas *c3 = new TCanvas("c3", "Conn.Comp. With Z", 300, 300 );
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
      color = selectCompColour(j, numColours);
      ConComps->SetMarkerColor(color);
      marker = (j % 34);
      if( marker < 20) { marker += 20; }
      if( marker > 34) { marker = 34; }
      ConComps->SetMarkerStyle(marker);
      ConComps->SetMarkerSize(0.8);
      // Conditions
      ePlotCondition  = "(EvtNum == "  + EeventNumString.str() + ") && ";
      ePlotCondition += "(CompNum == " + ecomponentNumber.str() + ")";
      ConComps->Draw("y:z_Det", ePlotCondition.c_str(), "same");
    }
    c3->Update();
    }
   
  // MC data plots.
     {
    TCanvas *c4 = new TCanvas("c4", "MC-Plot", 300, 300 );
    grid->Draw("y:x","","");
    TH2F *htemp = (TH2F*)gPad->GetPrimitive("htemp");
    htemp->GetXaxis()->SetTitle("x [cm]");
    htemp->GetYaxis()->SetTitle("y [cm]");
    htemp->SetTitle("");
    pos->SetMarkerSize(0.5);
    pos->Draw("my:mx","my > -10000","same");
    c4->Update();
  }
  // Read data
     {
    TCanvas *c5 = new TCanvas("c5", "Read-Plot", 300, 300 );
    grid->Draw("y:x","","");
    TH2F *htemp = (TH2F*)gPad->GetPrimitive("htemp");
    htemp->GetXaxis()->SetTitle("x [cm]");
    htemp->GetYaxis()->SetTitle("y [cm]");
    htemp->SetTitle("");
    pos->SetMarkerSize(0.5);
    pos->Draw("y:x","my > -10000","same");
    c5->Update();
  }
  // MC in Z-coordinate
      {
    TCanvas *c6 = new TCanvas("c6", "MC_Z_Plot", 300, 300 );
    grid->Draw("y:z","","");
    TH2F *htemp = (TH2F*)gPad->GetPrimitive("htemp");
    htemp->GetXaxis()->SetTitle("z [cm]");
    htemp->GetYaxis()->SetTitle("y [cm]");
    htemp->SetTitle("");
    pos->SetMarkerSize(0.5);
    pos->Draw("my:mz","my > -10000","same");
    c6->Update();
    }
    // MC in Z-coordinate
  {
    TCanvas *c8 = new TCanvas("c8", "Read_Z_Plot", 300, 300 );
    grid->Draw("y:z","","");
    TH2F *htemp = (TH2F*)gPad->GetPrimitive("htemp");
    htemp->GetXaxis()->SetTitle("z [cm]");
    htemp->GetYaxis()->SetTitle("y [cm]");
    htemp->SetTitle("");
    pos->SetMarkerSize(0.5);
    pos->Draw("y:z","my > -10000","same");
    c8->Update();
    }
    #endif
  
  //_____________ END TEMPORARY PLOTS May be deleted.
  c1->Divide(2,2);
  // Read data
  c1->cd(1);
  if(dim == 2) {
    grid->Draw("y:x","","");
    TH2F *htemp = (TH2F*)gPad->GetPrimitive("htemp");
    htemp->GetXaxis()->SetTitle("x [cm]");
    htemp->GetYaxis()->SetTitle("y [cm]");
    htemp->SetTitle("");
    pos->Draw("y:x","my > -10000","same");
  }
  if(dim == 3) {
    grid->Draw("y:x:z","","");
    TH3F *htemp = (TH3F*)gPad->GetPrimitive("htemp");
    htemp->GetXaxis()->SetTitle("z [cm]");
    htemp->GetYaxis()->SetTitle("x [cm]");
    htemp->GetZaxis()->SetTitle("y [cm]");
    htemp->SetTitle("");
    pos->Draw("y:x:z","my > -10000","same");
  }
  // MC data
  c1->cd(2);
  if(dim == 2) {
    grid->Draw("y:x","","");
    TH2F *htemp = (TH2F*)gPad->GetPrimitive("htemp");
    htemp->GetXaxis()->SetTitle("x [cm]");
    htemp->GetYaxis()->SetTitle("y [cm]");
    htemp->SetTitle("");
    pos->Draw("my:mx","my > -10000","same");
  }
  if(dim == 3) {
    grid->Draw("y:x:z","","");
    TH3F *htemp = (TH3F*)gPad->GetPrimitive("htemp");
    htemp->GetXaxis()->SetTitle("z [cm]");
    htemp->GetYaxis()->SetTitle("x [cm]");
    htemp->GetZaxis()->SetTitle("y [cm]");
    htemp->SetTitle("");
    pos->Draw("my:mx:mz","my > -10000","same");
  }
  
  // Connected components
  c1->cd(3);
  if(dim == 2) {
    grid->Draw("y:x","","");
    TH2F *htemp = (TH2F*)gPad->GetPrimitive("htemp");
    htemp->GetXaxis()->SetTitle("x [cm]");
    htemp->GetYaxis()->SetTitle("y [cm]");
    htemp->SetTitle("");
  }
  else {
    grid->Draw("y:x:z","","");
    TH3F *htemp = (TH3F*)gPad->GetPrimitive("htemp");
    htemp->GetXaxis()->SetTitle("z [cm]");
    htemp->GetYaxis()->SetTitle("x [cm]");
    htemp->GetZaxis()->SetTitle("y [cm]");
    htemp->SetTitle("");
  }
  
  std::stringstream EventNumString;
  EventNumString << evtNum;
  std::string PlotCondition;
  for(size_t j = 0; j < nComp; j++) {
    std::stringstream componentNumber;
    componentNumber << j;
    // Set marker and colours per component.
    color = selectCompColour(j, numColours);
    ConComps->SetMarkerColor(color);
    marker = (j % 34);
    if( marker < 20) { marker += 20; }
    if( marker > 34) { marker = 34; }
    ConComps->SetMarkerStyle(marker);
    ConComps->SetMarkerSize(0.8);
    // Conditions
    PlotCondition  = "(EvtNum == "  + EventNumString.str() + ") && ";
    PlotCondition += "(CompNum == " + componentNumber.str() + ")";
    if(dim == 2) {
      ConComps->Draw("y:x", PlotCondition.c_str(), "same");
    }
    if(dim == 3){
      ConComps->Draw("y:x:z", PlotCondition.c_str(), "same");
    }
  }

  // Merged and Z determined
  c1->cd(4);
  if(dim == 2) {
    grid->Draw("y:x","","");
    TH2F *htemp = (TH2F*)gPad->GetPrimitive("htemp");
    htemp->GetXaxis()->SetTitle("x [cm]");
    htemp->GetYaxis()->SetTitle("y [cm]");
    htemp->SetTitle("");
  }
  else {
    grid->Draw("y:x:z","","");
    TH2F *htemp = (TH2F*)gPad->GetPrimitive("htemp");
    htemp->GetXaxis()->SetTitle("z [cm]");
    htemp->GetYaxis()->SetTitle("x [cm]");
    htemp->GetZaxis()->SetTitle("y [cm]");
    htemp->SetTitle("");
  }
   
  for(size_t k = 0; k < nComp; k++) {
    std::stringstream componentNumber;
    componentNumber << k;
    // Select marker and color per component
    color = selectCompColour(k, numColours);
    ConComps->SetMarkerColor(color);
    marker = k % 34;    
    if( marker < 20) { marker += 20; }
    if( marker > 34) { marker = 34; }
    ConComps->SetMarkerStyle( marker );
    ConComps->SetMarkerSize(0.8);
    // Conditions
    PlotCondition  = "(EvtNum == "  + EventNumString.str() + ") && ";
    PlotCondition += "(CompNum == " + componentNumber.str() + ")";
    if(dim == 2) {
      ConComps->Draw("y:x", PlotCondition.c_str(), "same");
    }
    if(dim == 3){
      ConComps->Draw("y:x:z_Det", PlotCondition.c_str(), "same");
    }
  }
  // Save output to a file
  // TDirectory* direct = 0;
  /* direct = inp.GetDirectory(dirName.c_str());
  if( !direct) {
    inp.mkdir(dirName.c_str());
  }
  inp.cd(dirName.c_str());
  c1->Write();*/
}
