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
#include "TMarker.h"

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
				  Double_t const ww = 600.00,
				  Double_t const hh = 600.00,
				  size_t  const dim = 2,
				  size_t evtNum = 0)
{
  TFile inp(InoutFile.c_str(),"UPDATE");

  // Get COnnected components
  TNtuple* ConComps = (TNtuple*) inp.Get("ConnectedComponents");
  TNtuple* ConCompsAnc = (TNtuple*) inp.Get("ConnectedComponentsAnchors");

  TNtuple* Poss   = (TNtuple*) inp.Get("Pos");
  TNtuple* MCposs = (TNtuple*) inp.Get("MCpos");

  // Polar coordinates.
  TNtuple* MCposPolar = (TNtuple*) inp.Get("MCposPolar");
  TNtuple* PosPolar   = (TNtuple*) inp.Get("PosPolar");
  TNtuple *ErrorPT = (TNtuple*) inp.Get("PerTrackError");

  Poss->AddFriend("MCpos");

  PosPolar->AddFriend("MCposPolar");

  // Get number of components for the current event
  float numComponents;

  TTree* EvtNumComps = (TTree*) inp.Get("ComponentPerEvt");
  EvtNumComps->SetBranchAddress("numComponents", &numComponents);

  EvtNumComps->GetEntry(evtNum);
  unsigned int nComp = static_cast<unsigned int>(numComponents);

  // Fetch data and make plots
  TNtuple* grid = (TNtuple*) inp.Get("ExtendedGrid");
  int numColours = 8;
  int color = 1;
  grid->SetMarkerSize(0.3);
  // Change to directory with input events
  inp.cd("InputEvents");

  std::stringstream evN;
  std::string tupName = "Evt";
  evN << evtNum;
  tupName = tupName + "_" + evN.str() + "_CoordsTuple";
  tupName = "InputEvents/" + tupName;

  // float *trackID;
  float_t id;
  std::vector< int > idtracks;
 
  TTree* tracknum = (TTree*) inp.Get(tupName.c_str());
  TBranch *bpx = tracknum->GetBranch("trackID");
  bpx->SetAddress(&id);
  
  Int_t nevent = (Int_t)tracknum->GetEntries();
  for (Int_t i=0;i<nevent;i++) {    
    bpx->GetEntry(i); //read branch fTracks.fPx
 
    // printf("cur track=%f \n",id);
    int intid =  static_cast<unsigned int>(id);
    if(!(std::find(idtracks.begin(), idtracks.end(), intid) != idtracks.end()))
      	idtracks.push_back(intid);

  }

  float_t idt;
  float_t comp;

  std::vector< int > idMatch;

  //TTree* tracknum = (TTree*) inp.Get(tupName.c_str());
  TBranch *bpxx = ConComps->GetBranch("bestIdx");
  TBranch *bpxxx = ConComps->GetBranch("CompNum");

  bpxx->SetAddress(&idt);
  bpxxx->SetAddress(&comp);

  int lastcc = -1;
  nevent = (Int_t)ConComps->GetEntries();
  for (Int_t i=0;i<nevent;i++) {    
    bpxx->GetEntry(i); //read branch fTracks.fPx
    bpxxx->GetEntry(i); //read branch fTracks.fPx

    // printf("ntracks=%f \n",id);
    int curcomp = static_cast<unsigned int>(comp);
    // printf("Comp %d\n",curcomp);

    if(curcomp != lastcc){
      int intidm =  static_cast<unsigned int>(idt);
      //  printf("%d\n",intidm);
      idMatch.push_back(intidm);

      lastcc = curcomp;
    }

     //  if(!(std::find(idtracks.begin(), idtracks.end(), intid) != idtracks.end()))
      //	idtracks.push_back(intid);
  }

  //printf("%d\n", idtracks.size());

  
  TNtuple* pos = (TNtuple*) inp.Get(tupName.c_str());

  
  // Create and set the canvas name and properties
  std::string canvasName = "Event_" + evN.str() + "_PlotCanvas";
  std::string descr = "All Collected Coordinates Event " + evN.str();

  // TCanvas *c1 = new TCanvas(canvasName.c_str(), descr.c_str(), ww, hh );

  int marker = 0;
  int badmarker = 24;
  float sizelabels = 0.05;
  float padleft = 0.1;
  float padtop = 0.05;
  float markersize = 0.8;
  //_____________ BEGIN TEMPORARY PLOTS May be DELETED
#if(CREATE_SEPARATE_PLOTS > 0)
  {
    // Connected components
    TCanvas *c2 = new TCanvas("c2", "Connected Components", ww, hh );
    c2->SetRightMargin(padtop);
    c2->SetLeftMargin(padleft);
    c2->SetTopMargin(padtop);
    c2->SetBottomMargin(padleft);
    grid->Draw("y:x","","");
    TH2F *htemp = (TH2F*)gPad->GetPrimitive("htemp");
    htemp->GetXaxis()->SetTitle("x [cm]");
    htemp->GetYaxis()->SetTitle("y [cm]");
    htemp->GetXaxis()->SetLabelSize(sizelabels);
    htemp->GetYaxis()->SetLabelSize(sizelabels);
    htemp->GetXaxis()->SetTitleSize(sizelabels);
    htemp->GetYaxis()->SetTitleOffset(1);
    htemp->GetYaxis()->SetTitleSize(sizelabels);

    htemp->SetTitle("");
    std::stringstream EeventNumString;
    EeventNumString << evtNum;
    std::string ePlotCondition;
    for(size_t j = 0; j < nComp; j++) {
      std::stringstream ecomponentNumber;
      ecomponentNumber << nComp-1-j;
      // Set marker and colours per component.
      int idmctr = idMatch[ nComp-1-j];
      marker = 20;//(j % 34);

      if(idmctr != -1)
	color = selectCompColour(idmctr, numColours);
      else{
	color = 1;
	marker = badmarker++;
      }
      ConComps->SetMarkerColor(color);
      //if( marker < 20) { marker += 20; }
      //if( marker > 34) { marker = 34; }
      ConComps->SetMarkerStyle(marker);
      ConComps->SetMarkerSize(markersize);
      // Conditions
      ePlotCondition  = "(EvtNum == "  + EeventNumString.str() + ") && ";
      ePlotCondition += "(CompNum == " + ecomponentNumber.str() + ")";
      ConComps->Draw("y_Det:x_Det", ePlotCondition.c_str(), "same");
      /*TMarker *m = new TMarker(24, 20, 29);
      m->SetMarkerSize(2.);
      m->Draw();
      TMarker *m1 = new TMarker(-10.5, -21, 29);
      m1->SetMarkerSize(2.);
      m1->Draw();*/
    }
    c2->Update();
    c2->SaveAs("myalgo_cm_xy.pdf");

    }
  // Z-reconstructed Connected components
   {
    TCanvas *c3 = new TCanvas("c3", "Conn.Comp. With Z", ww, hh );
    c3->SetRightMargin(padtop);
    c3->SetLeftMargin(padleft);
    c3->SetTopMargin(padtop);
    c3->SetBottomMargin(padleft);
    grid->Draw("y:z","","");
    TH2F *htemp = (TH2F*)gPad->GetPrimitive("htemp");
    
    htemp->GetXaxis()->SetTitle("z [cm]");
    htemp->GetYaxis()->SetTitle("y [cm]");
        htemp->GetXaxis()->SetLabelSize(sizelabels);
    htemp->GetYaxis()->SetLabelSize(sizelabels);
    htemp->GetXaxis()->SetTitleSize(sizelabels);
    htemp->GetYaxis()->SetTitleOffset(1);
    htemp->GetYaxis()->SetTitleSize(sizelabels);
    htemp->SetTitle("");
    std::stringstream EeventNumString;
    EeventNumString << evtNum;
    std::string ePlotCondition;
    for(size_t j = 0; j < nComp; j++) {
      std::stringstream ecomponentNumber;
      
      ecomponentNumber << j;
      std::stringstream edetID;
      edetID << 5000; 
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
      ConComps->SetMarkerSize(markersize);
      // Conditions
         ePlotCondition  = "(EvtNum == "  + EeventNumString.str() + ") && ";
      ePlotCondition += "(CompNum == " + ecomponentNumber.str() + ") ";
      //ePlotCondition += "(z_Det !=  0 ) &&";
      //ePlotCondition += "(tubeId < " + edetID.str() + ")";

      ConComps->Draw("y_Det:z_Det", ePlotCondition.c_str(), "same");
  

      /*  ConComps->SetMarkerStyle(29);
      ConComps->SetMarkerSize(2.);
      ePlotCondition  = "(EvtNum == "  + EeventNumString.str() + ") && ";
      ePlotCondition += "(CompNum == " + ecomponentNumber.str() + ") && ";
      ePlotCondition += "(tubeId >= " + edetID.str() + ")";
      ConComps->Draw("y_Det:z_Det", ePlotCondition.c_str(), "same");*/
    }
    c3->Update();
    c3->SaveAs("myalgo_cm_z.pdf");

    }
   
  // MC data plots.
  {
    TCanvas *c4 = new TCanvas("c4", "MC-Plot", ww, hh );
    c4->SetRightMargin(padtop);
    c4->SetLeftMargin(padleft);
    c4->SetTopMargin(padtop);
    c4->SetBottomMargin(padleft);
    grid->Draw("y:x","","");
    TH2F *htemp = (TH2F*)gPad->GetPrimitive("htemp");
    htemp->GetXaxis()->SetTitle("x [cm]");
    htemp->GetYaxis()->SetTitle("y [cm]");
        htemp->GetXaxis()->SetLabelSize(sizelabels);
    htemp->GetYaxis()->SetLabelSize(sizelabels);
    htemp->GetXaxis()->SetTitleSize(sizelabels);
    htemp->GetYaxis()->SetTitleOffset(1);
    htemp->GetYaxis()->SetTitleSize(sizelabels);
    htemp->SetTitle("");
    std::stringstream EeventNumString;
    EeventNumString << evtNum;
    std::string ePlotCondition;
    for(size_t j = 0; j < idtracks.size(); j++) {
      std::stringstream ecomponentNumber;
      ecomponentNumber << idtracks[j];
      int idmctr = idtracks[j];
      //printf("%d\n", idmctr);

      color = selectCompColour(idmctr, numColours);
      // printf("%d, %d\n", j, color);

      pos->SetMarkerColor(color);
      marker = 20;//(j % 34);
      // if( marker < 20) { marker += 20; }
      //if( marker > 34) { marker = 34; }
      pos->SetMarkerStyle(marker);
      pos->SetMarkerSize(markersize);
      ePlotCondition  = "(EvtNum == "  + EeventNumString.str() + ") && ";
      ePlotCondition += "(trackID == " + ecomponentNumber.str() + ") && my > -10000";
      pos->Draw("my:mx",ePlotCondition.c_str(),"same");
    }
    c4->Update();
     c4->SaveAs("MC_xy.pdf");

  }
  // Read data
   {
    TCanvas *c5 = new TCanvas("c5", "Read-Plot", ww, hh );
    c5->SetRightMargin(padtop);
    c5->SetLeftMargin(padleft);
    c5->SetTopMargin(padtop);
    c5->SetBottomMargin(padleft);
    grid->Draw("y:x","","");
    TH2F *htemp = (TH2F*)gPad->GetPrimitive("htemp");
    htemp->GetXaxis()->SetTitle("x [cm]");
    htemp->GetYaxis()->SetTitle("y [cm]");
        htemp->GetXaxis()->SetLabelSize(sizelabels);
    htemp->GetYaxis()->SetLabelSize(sizelabels);
    htemp->GetXaxis()->SetTitleSize(sizelabels);
    htemp->GetYaxis()->SetTitleOffset(1);
    htemp->GetYaxis()->SetTitleSize(sizelabels);
    htemp->SetTitle("");
    std::stringstream EeventNumString;
    EeventNumString << evtNum;
    std::string ePlotCondition;
    for(size_t j = 0; j < idtracks.size(); j++) {
      std::stringstream ecomponentNumber;
      ecomponentNumber << idtracks[j];
      color = selectCompColour(j, numColours);
      pos->SetMarkerColor(color);
      marker = 20;//(j % 34);
      // if( marker < 20) { marker += 20; }
      //if( marker > 34) { marker = 34; }
      pos->SetMarkerStyle(marker);
      pos->SetMarkerSize(markersize);
      ePlotCondition  = "(EvtNum == "  + EeventNumString.str() + ") && ";
      ePlotCondition += "(trackID == " + ecomponentNumber.str() + ") && my > -10000";
      pos->Draw("y:x",ePlotCondition.c_str(),"same");
    }
    c5->Update();
     c5->SaveAs("Read_xy.pdf");

  }
  // MC in Z-coordinate
  {
    TCanvas *c6 = new TCanvas("c6", "MC_Z_Plot", ww, hh );
    c6->SetRightMargin(padtop);
    c6->SetLeftMargin(padleft);
    c6->SetTopMargin(padtop);
    c6->SetBottomMargin(padleft);
    grid->Draw("y:z","","");
    TH2F *htemp = (TH2F*)gPad->GetPrimitive("htemp");
    htemp->GetXaxis()->SetTitle("z [cm]");
    htemp->GetYaxis()->SetTitle("y [cm]");
        htemp->GetXaxis()->SetLabelSize(sizelabels);
    htemp->GetYaxis()->SetLabelSize(sizelabels);
    htemp->GetXaxis()->SetTitleSize(sizelabels);
    htemp->GetYaxis()->SetTitleOffset(1);
    htemp->GetYaxis()->SetTitleSize(sizelabels);
    htemp->SetTitle("");
    std::stringstream EeventNumString;
    EeventNumString << evtNum;
    std::string ePlotCondition;
    for(size_t j = 0; j < idtracks.size(); j++) {
      std::stringstream ecomponentNumber;
      ecomponentNumber << idtracks[j];
      color = selectCompColour(j, numColours);
      pos->SetMarkerColor(color);
      marker = 20;//(j % 34);
      // if( marker < 20) { marker += 20; }
      //if( marker > 34) { marker = 34; }
      pos->SetMarkerStyle(marker);
      pos->SetMarkerSize(0.8);
      ePlotCondition  = "(EvtNum == "  + EeventNumString.str() + ") && ";
      ePlotCondition += "(trackID == " + ecomponentNumber.str() + ") && my > -10000";
      pos->Draw("my:mz",ePlotCondition.c_str(),"same");
    }
    c6->Update();
    c6->SaveAs("MC_yz.pdf");

  }

  if(dim == 3) {
    TCanvas *c11 = new TCanvas("c11", "MC_Plot3D", ww, hh );
    grid->Draw("y:x:z","","");
    TH3F *htemp = (TH3F*)gPad->GetPrimitive("htemp");
    htemp->GetXaxis()->SetTitle("z [cm]");
    htemp->GetYaxis()->SetTitle("x [cm]");
    htemp->GetZaxis()->SetTitle("y [cm]");
    htemp->GetXaxis()->SetLabelSize(0.05);
    htemp->GetYaxis()->SetLabelSize(0.05);
    htemp->GetXaxis()->SetTitleSize(.05);
    htemp->GetYaxis()->SetTitleOffset(1);
    htemp->GetYaxis()->SetTitleSize(.05);
    htemp->SetTitle("");
    std::stringstream EeventNumString;
    EeventNumString << evtNum;
    std::string ePlotCondition;
    for(size_t j = 0; j < idtracks.size(); j++) {
      std::stringstream ecomponentNumber;
      ecomponentNumber << idtracks[j];
      color = selectCompColour(j, numColours);
      pos->SetMarkerColor(color);
      marker = 20;//(j % 34);
      // if( marker < 20) { marker += 20; }
      //if( marker > 34) { marker = 34; }
      pos->SetMarkerStyle(marker);
      pos->SetMarkerSize(0.8);
      ePlotCondition  = "(EvtNum == "  + EeventNumString.str() + ") && ";
      ePlotCondition += "(trackID == " + ecomponentNumber.str() + ") && my > -10000";
      pos->Draw("my:mx:mz",ePlotCondition.c_str(),"same");
    }

    if(dim == 3) {
      TCanvas *c12 = new TCanvas("c12", "Det_Plot3D", 300, 300 );
      grid->Draw("y:x:z","","");
      TH3F *htemp = (TH3F*)gPad->GetPrimitive("htemp");
      htemp->GetXaxis()->SetTitle("z [cm]");
      htemp->GetYaxis()->SetTitle("x [cm]");
      htemp->GetZaxis()->SetTitle("y [cm]");
      htemp->SetTitle("");
      std::stringstream EeventNumString;
      EeventNumString << evtNum;
      std::string ePlotCondition;
      for(size_t j = 0; j < idtracks.size(); j++) {
	std::stringstream ecomponentNumber;
	ecomponentNumber << idtracks[j];
	color = selectCompColour(j, numColours);
	pos->SetMarkerColor(color);
	marker = 20;//(j % 34);
      // if( marker < 20) { marker += 20; }
      //if( marker > 34) { marker = 34; }
	pos->SetMarkerStyle(marker);
	pos->SetMarkerSize(0.8);
	ePlotCondition  = "(EvtNum == "  + EeventNumString.str() + ") && ";
	ePlotCondition += "(trackID == " + ecomponentNumber.str() + ") && my > -10000";
	pos->Draw("y:x:z",ePlotCondition.c_str(),"same");
      }
    }

    
    if(dim == 3) {
      TCanvas *c13 = new TCanvas("c13", "CC_Plot3D", 300, 300 );
      grid->Draw("y:x:z","","");
      TH3F *htemp = (TH3F*)gPad->GetPrimitive("htemp");
      htemp->GetXaxis()->SetTitle("z [cm]");
      htemp->GetYaxis()->SetTitle("x [cm]");
      htemp->GetZaxis()->SetTitle("y [cm]");
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
	if(idmctr != -1)
	  color = selectCompColour(idmctr, numColours);
	else{
	  color = 1;
	  marker = badmarker++;;
	}
      //	color = selectCompColour(j, numColours);
	ConComps->SetMarkerColor(color);
	//marker = 20;//(j % 34);
      // if( marker < 20) { marker += 20; }
      //if( marker > 34) { marker = 34; }
	ConComps->SetMarkerStyle(marker);
	ConComps->SetMarkerSize(0.8);
	// Conditions
	ePlotCondition  = "(EvtNum == "  + EeventNumString.str() + ") && ";
	ePlotCondition += "(CompNum == " + ecomponentNumber.str() + ")";
	ConComps->Draw("y_Det:x_Det:z_Det",ePlotCondition.c_str(),"same");
      }
    }
    
  }
	  
  

  {


    // Connected components
    TCanvas *c20 = new TCanvas("c20", "Connected Component Anchors", 300, 300 );
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
      ConCompsAnc->SetMarkerColor(color);
      marker = (j % 34);
      if( marker < 20) { marker += 20; }
      if( marker > 34) { marker = 34; }
      ConCompsAnc->SetMarkerStyle(marker);
      ConCompsAnc->SetMarkerSize(0.8);
      // Conditions
      ePlotCondition  = "(EvtNum == "  + EeventNumString.str() + ") && ";
      ePlotCondition += "(CompNum == " + ecomponentNumber.str() + ")";
      ConCompsAnc->Draw("y_Det:x_Det", ePlotCondition.c_str(), "same");
    }
    c20->Update();
    // c2->SaveAs("myalgo_cm_xy.pdf");

  }/*
  // Z-reconstructed Connected components
  {
    TCanvas *c30 = new TCanvas("c30", "Conn.Comp. With Z2", 300, 300 );
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
      ConCompsAnc->SetMarkerColor(color);
      marker = (j % 34);
      if( marker < 20) { marker += 20; }
      if( marker > 34) { marker = 34; }
      ConCompsAnc->SetMarkerStyle(marker);
      ConCompsAnc->SetMarkerSize(0.8);
      // Conditions
      ePlotCondition  = "(EvtNum == "  + EeventNumString.str() + ") && ";
      ePlotCondition += "(CompNum == " + ecomponentNumber.str() + ")";
      ConCompsAnc->Draw("y_Det:z_Det", ePlotCondition.c_str(), "same");
    }
    c30->Update();
    //c3->SaveAs("myalgo_cm_z.pdf");

    }*/
  
    #endif



  /*
  
  //_____________ END TEMPORARY PLOTS May be deleted.
  c1->Divide(3,3);
  // Read data
  c1->cd(1);
  if(dim >= 2) {
    grid->Draw("y:x","","");
    TH2F *htemp = (TH2F*)gPad->GetPrimitive("htemp");
    htemp->GetXaxis()->SetTitle("x [cm]");
    htemp->GetYaxis()->SetTitle("y [cm]");
    htemp->SetTitle("");
        std::stringstream EeventNumString;
    EeventNumString << evtNum;
    std::string ePlotCondition;
    for(size_t j = 0; j < idtracks.size(); j++) {
      std::stringstream ecomponentNumber;
      ecomponentNumber << idtracks[j];
      color = selectCompColour(j, numColours);
      pos->SetMarkerColor(color);
      marker = 20;
      pos->SetMarkerStyle(marker);
      pos->SetMarkerSize(0.8);
      ePlotCondition  = "(EvtNum == "  + EeventNumString.str() + ") && ";
      ePlotCondition += "(trackID == " + ecomponentNumber.str() + ") && my > -10000";
      pos->Draw("my:mx",ePlotCondition.c_str(),"same");
    }
  } 
  c1->cd(2);

  if(dim >= 2) {
    grid->Draw("y:z","","");
    TH2F *htemp = (TH2F*)gPad->GetPrimitive("htemp");
    htemp->GetXaxis()->SetTitle("z [cm]");
    htemp->GetYaxis()->SetTitle("y [cm]");
    htemp->SetTitle("");
        std::stringstream EeventNumString;
    EeventNumString << evtNum;
    std::string ePlotCondition;
    for(size_t j = 0; j < idtracks.size(); j++) {
      std::stringstream ecomponentNumber;
      ecomponentNumber << idtracks[j];
      color = selectCompColour(j, numColours);
      pos->SetMarkerColor(color);
      marker = 20;
      pos->SetMarkerStyle(marker);
      pos->SetMarkerSize(0.8);
      ePlotCondition  = "(EvtNum == "  + EeventNumString.str() + ") && ";
      ePlotCondition += "(trackID == " + ecomponentNumber.str() + ") && my > -10000";
      pos->Draw("my:mz",ePlotCondition.c_str(),"same");
    }
   }
   
  // MC data

  c1->cd(3);
  c1->DrawFrame(15,-190,42,190);
 if(dim >= 2) {
    grid->Draw("","","");
    TH2F *htemp = (TH2F*)gPad->GetPrimitive("htemp");
    //htemp->Draw();
    // htemp->GetXaxis()->SetTitle("x [cm]");
    //htemp->GetYaxis()->SetTitle("y [cm]");
    // htemp->SetTitle("");
    std::stringstream EeventNumString;
    EeventNumString << evtNum;
    std::string ePlotCondition;
    for(size_t j = 0; j < idtracks.size(); j++) {
      std::stringstream ecomponentNumber;
      ecomponentNumber << idtracks[j];
      color = selectCompColour(j, numColours);
      pos->SetMarkerColor(color);
      marker = 20;
      pos->SetMarkerStyle(marker);
      pos->SetMarkerSize(0.8);
      ePlotCondition  = "(EvtNum == "  + EeventNumString.str() + ") && ";
      ePlotCondition += "(trackID == " + ecomponentNumber.str() + ") && my > -10000";
      // if (j == 0) pos->Draw("mthetaDeg:mr",ePlotCondition.c_str(),"");
      pos->Draw("mthetaDeg:mr",ePlotCondition.c_str(),"same");
    }
    }


  
  c1->cd(4);
  if(dim >= 2) {
    grid->Draw("y:x","","");
    TH2F *htemp = (TH2F*)gPad->GetPrimitive("htemp");
    htemp->GetXaxis()->SetTitle("x [cm]");
    htemp->GetYaxis()->SetTitle("y [cm]");
    htemp->SetTitle("");
        std::stringstream EeventNumString;
    EeventNumString << evtNum;
    std::string ePlotCondition;
    for(size_t j = 0; j < idtracks.size(); j++) {
      std::stringstream ecomponentNumber;
      ecomponentNumber << idtracks[j];
      color = selectCompColour(j, numColours);
      pos->SetMarkerColor(color);
      marker =20;
      pos->SetMarkerStyle(marker);
      pos->SetMarkerSize(0.8);
      ePlotCondition  = "(EvtNum == "  + EeventNumString.str() + ") && ";
      ePlotCondition += "(trackID == " + ecomponentNumber.str() + ") && my > -10000";
      pos->Draw("y:x",ePlotCondition.c_str(),"same");
    }
     }

  c1->cd(5);
   if(dim >= 2) {
    grid->Draw("y:z","","");
    TH2F *htemp = (TH2F*)gPad->GetPrimitive("htemp");
    htemp->GetXaxis()->SetTitle("z [cm]");
    htemp->GetYaxis()->SetTitle("y [cm]");
    htemp->SetTitle("");
        std::stringstream EeventNumString;
    EeventNumString << evtNum;
    std::string ePlotCondition;
    for(size_t j = 0; j < idtracks.size(); j++) {
      std::stringstream ecomponentNumber;
      ecomponentNumber << idtracks[j];
      color = selectCompColour(j, numColours);
      pos->SetMarkerColor(color);
      marker =20;
      pos->SetMarkerStyle(marker);
      pos->SetMarkerSize(0.8);
      ePlotCondition  = "(EvtNum == "  + EeventNumString.str() + ") && ";
      ePlotCondition += "(trackID == " + ecomponentNumber.str() + ") && my > -10000";
      pos->Draw("y:z",ePlotCondition.c_str(),"same");
    }
   }

   c1->cd(6);
  c1->DrawFrame(15,-190,42,190);
  if(dim >= 2) {
    grid->Draw("","","");
    TH2F *htemp = (TH2F*)gPad->GetPrimitive("htemp");
    //htemp->Draw();
    // htemp->GetXaxis()->SetTitle("x [cm]");
    //htemp->GetYaxis()->SetTitle("y [cm]");
    // htemp->SetTitle("");
    std::stringstream EeventNumString;
    EeventNumString << evtNum;
    std::string ePlotCondition;
    for(size_t j = 0; j < idtracks.size(); j++) {
      std::stringstream ecomponentNumber;
      ecomponentNumber << idtracks[j];
      color = selectCompColour(j, numColours);
      pos->SetMarkerColor(color);
      marker =20;
      pos->SetMarkerStyle(marker);
      pos->SetMarkerSize(0.8);
      ePlotCondition  = "(EvtNum == "  + EeventNumString.str() + ") && ";
      ePlotCondition += "(trackID == " + ecomponentNumber.str() + ") && my > -10000";
      // if (j == 0) pos->Draw("mthetaDeg:mr",ePlotCondition.c_str(),"");
      pos->Draw("thetaDeg:r",ePlotCondition.c_str(),"same");
    }
   }
 
  // Connected components
  c1->cd(7);
 if(dim >= 2) {
    grid->Draw("y:x","","");
    TH2F *htemp = (TH2F*)gPad->GetPrimitive("htemp");
    htemp->GetXaxis()->SetTitle("x [cm]");
    htemp->GetYaxis()->SetTitle("y [cm]");
    htemp->SetTitle("");
    
   }

  std::stringstream EventNumString;
  EventNumString << evtNum;
  std::string PlotCondition;
  for(size_t j = 0; j < nComp; j++) {
    std::stringstream componentNumber;
    componentNumber << j;
    // Set marker and colours per component.
     marker = 20;//(j % 34);
    int idmctr = idMatch[j];

    if(idmctr != -1)
      color = selectCompColour(idmctr, numColours);
    else{
      color = 1;
      marker = badmarker++;;
    }
    ConComps->SetMarkerColor(color);
    //    marker = 20;
    ConComps->SetMarkerStyle(marker);
    ConComps->SetMarkerSize(0.8);
    // Conditions
    PlotCondition  = "(EvtNum == "  + EventNumString.str() + ") && ";
    PlotCondition += "(CompNum == " + componentNumber.str() + ")";
    //if(dim == 2) {
      ConComps->Draw("y:x", PlotCondition.c_str(), "same");
      //}
  }

  // Merged and Z determined
  c1->cd(8);
  if(dim >= 2) {
    grid->Draw("y:z","","");
    TH2F *htemp = (TH2F*)gPad->GetPrimitive("htemp");
    htemp->GetXaxis()->SetTitle("z [cm]");
    htemp->GetYaxis()->SetTitle("y [cm]");
    htemp->SetTitle("");
  }
  
  for(size_t k = 0; k < nComp; k++) {
    std::stringstream componentNumber;
    componentNumber << k;
    // Select marker and color per component
    marker = 20;//(j % 34);
    int idmctr = idMatch[k];

    if(idmctr != -1)
      color = selectCompColour(idmctr, numColours);
    else{
      color = 1;
      marker = badmarker++;;
    }
    color = selectCompColour(k, numColours);
    ConComps->SetMarkerColor(color);
    ConComps->SetMarkerStyle( marker );
    ConComps->SetMarkerSize(0.8);
    // Conditions
    PlotCondition  = "(EvtNum == "  + EventNumString.str() + ") && ";
    PlotCondition += "(CompNum == " + componentNumber.str() + ")";
    //  if(dim == 2) {
      ConComps->Draw("y:z_Det", PlotCondition.c_str(), "same");
      //  }
    // if(dim == 3){
    //   ConComps->Draw("y:x:z_Det", PlotCondition.c_str(), "same");
      //  }
  }

  c1->cd(9);
  c1->DrawFrame(15,-190,42,190);

  if(dim >= 2) {
    //  grid->Draw("y:z","","");
    TH2F *htemp = (TH2F*)gPad->GetPrimitive("htemp");
    //   htemp->GetXaxis()->SetTitle("y [cm]");
    // htemp->GetYaxis()->SetTitle("z [cm]");
    //  htemp->SetTitle("");
  }
   
  for(size_t k = 0; k < nComp; k++) {
    std::stringstream componentNumber;
    componentNumber << k;
    // Select marker and color per component
    color = selectCompColour(k, numColours);
    ConComps->SetMarkerColor(color);
    marker = 20;
    ConComps->SetMarkerStyle( marker );
    ConComps->SetMarkerSize(0.8);
    // Conditions
    PlotCondition  = "(EvtNum == "  + EventNumString.str() + ") && ";
    PlotCondition += "(CompNum == " + componentNumber.str() + ")";
    //  if(dim == 2) {
      ConComps->Draw("thetaDeg:r", PlotCondition.c_str(), "same");
      // }
  }
  // Save output to a file
  // TDirectory* direct = 0;
  direct = inp.GetDirectory(dirName.c_str());
  if( !direct) {
    inp.mkdir(dirName.c_str());
  }
  inp.cd(dirName.c_str());
  c1->Write();*/
}
