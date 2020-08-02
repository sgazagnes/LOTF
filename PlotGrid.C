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

#define SAVE_PDF_PLOT 1

void PlotGrid( std::string const& Infile = "Tracks_output.root",
	       Double_t const ww      = 300,
	       Double_t const hh      = 300,
	       size_t dim = 2
	       )
{
  TFile inp(Infile.c_str(),"READ");
  //TCanvas *GridPoints = new TCanvas("GridPoints", "GridPoints", ww, hh );
  
  TH2F *bla = new TH2F("Components","Components",80,-40,40, 80,-40,40);
  bla->Draw();
  
  // Original Grid
  TNtuple* plot = (TNtuple*) inp.Get("OrigGridCoord");
  plot->SetMarkerColor(14);
  plot->SetMarkerStyle(6);
  plot->SetMarkerSize(0.6);
  if(dim == 2) {
    plot->Draw("y:x","","");
  }
  else {
    plot->Draw("y:x:z","","");
  }

  // Virtual nodes
  TNtuple* vTube = (TNtuple*) inp.Get("VirtualNodes");
  vTube->SetMarkerColor(6);
  vTube->SetMarkerStyle(6);
  vTube->SetMarkerSize(0.6);
  if(dim == 2) {
    vTube->Draw("y:x","","same");
  }
  else {
    vTube->Draw("y:x:z","","same");
  }

  //_____________ BEGIN TEMPORARY PLOTS
  {
    TCanvas *c2 = new TCanvas("TmporaryMC", "blabla", 300, 300 );
    plot->Draw("y:x","","");
    vTube->Draw("y:x","","same");
    
    TH2F *htemp = (TH2F*)gPad->GetPrimitive("htemp");
    htemp->GetXaxis()->SetTitle("x [cm]");
    htemp->GetYaxis()->SetTitle("y [cm]");
    htemp->SetTitle("");
    c2->Update();
  }
  {
    TCanvas *c3 = new TCanvas("TmporaryREAD", "blabla", 300, 300 );
    plot->Draw("y:z","","");
    vTube->Draw("y:z","","same");
    
    TH2F *htemp = (TH2F*)gPad->GetPrimitive("htemp");
    htemp->GetXaxis()->SetTitle("z [cm]");
    htemp->GetYaxis()->SetTitle("y [cm]");
    htemp->SetTitle("");
    c3->Update();
  }
  //_____________ END TEMPORARY PLOTS
}
