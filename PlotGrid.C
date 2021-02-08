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
float theta =-90;
float phi = 180;
/*void Animate()
{

   theta += 0;
   phi   += 0.5;
   gPad->GetView()->RotateView(theta,phi);
   gPad->Modified();
   gPad->Update();
}

void Animate2()
{

   theta += 2;
   phi   += 0;
   gPad->GetView()->RotateView(theta,phi);
   gPad->Modified();
   gPad->Update();
   }*/
void PlotGrid( std::string const& Infile = "Tracks_output.root",
	       Double_t const ww      = 800,
	       Double_t const hh      = 800,
	       size_t dim = 3
	       )
{
  TFile inp(Infile.c_str(),"READ");
  TCanvas *GridPoints = new TCanvas("GridPoints", "GridPoints", ww, hh );
  
  /*TH2F *bla = new TH2F("Components","Components",80,-40,40, 80,-40,40);
  bla->GetXaxis()->SetTitle("x [cm]");
  bla->GetYaxis()->SetTitle("y [cm]");
  bla->SetTitle("");
  bla->Draw();*/
  /* TH2F *htemp = (TH2F*)gPad->GetPrimitive("htemp");
  htemp->GetXaxis()->SetTitle("x [cm]");
  htemp->GetYaxis()->SetTitle("y [cm]");
  htemp->SetTitle("");
  */  // Original Grid
  TNtuple* plot = (TNtuple*) inp.Get("OrigGridCoord");
  plot->SetMarkerColor(kBlack);
  plot->SetMarkerStyle(6);
  plot->SetMarkerSize(1);
  if(dim == 2) {
    plot->Draw("y:x","",""); 
  }
  else {
    plot->Draw("y:x:z","","");
  }

  // Virtual nodes
  TNtuple* vTube = (TNtuple*) inp.Get("VirtualNodesLayer");
  vTube->SetMarkerColor(kMagenta);
  vTube->SetMarkerStyle(6);
  vTube->SetMarkerSize(1);
  if(dim == 2) {
    vTube->Draw("y:x","","same");
  }
  else {
    vTube->Draw("y:x:z","","same");
  }

  TNtuple* vTubeS = (TNtuple*) inp.Get("VirtualNodesSector");
  vTubeS->SetMarkerColor(8);
  vTubeS->SetMarkerStyle(6);
  vTubeS->SetMarkerSize(1);
  if(dim == 2) {
    vTubeS->Draw("y:x","","same");
  }
  else {
    vTubeS->Draw("y:x:z","","same");
  }
  TH2F *htemp = (TH2F*)gPad->GetPrimitive("htemp");
  htemp->GetXaxis()->SetTitle("z [cm]");
  htemp->GetYaxis()->SetTitle("x [cm]");
  htemp->GetZaxis()->SetTitle("y [cm]");

  htemp->SetTitle("");
  GridPoints->Update();
  GridPoints->SaveAs("grid_xyz.pdf");

  
  //_____________ BEGIN TEMPORARY PLOTS
   {
    TCanvas *c2 = new TCanvas("TmporaryMC", "blabla",  ww, hh );
    plot->Draw("y:x","","");
    // vTube->Draw("y:x","","same");
    // vTubeS->Draw("y:x","","same");

    TH2F *htemp = (TH2F*)gPad->GetPrimitive("htemp");
    htemp->GetXaxis()->SetTitle("x [cm]");
    htemp->GetYaxis()->SetTitle("y [cm]");
    htemp->SetTitle("");
    c2->Update();
    // c2->SaveAs("grid_xy.pdf");
  }
  {
    TCanvas *c3 = new TCanvas("TmporaryREAD", "blabla",  ww, hh );
    plot->Draw("y:z","","");
    //vTube->Draw("y:z","","same");
    // vTubeS->Draw("y:z","","same");

    TH2F *htemp = (TH2F*)gPad->GetPrimitive("htemp");
    htemp->GetXaxis()->SetTitle("z [cm]");
    htemp->GetYaxis()->SetTitle("y [cm]");
    htemp->SetTitle("");
    c3->Update();
    //  c3->SaveAs("grid_yz.pdf");

    }

  {
    TCanvas *c4 = new TCanvas("TmporaryREAD2", "blabla", ww, hh );
    plot->Draw("x:z","","");
    //    vTube->Draw("x:z","","same");
    //  vTubeS->Draw("x:z","","same");

    TH2F *htemp = (TH2F*)gPad->GetPrimitive("htemp");
    htemp->GetXaxis()->SetTitle("z [cm]");
    htemp->GetYaxis()->SetTitle("x [cm]");
    htemp->SetTitle("");
    c4->Update();
    // c4->SaveAs("grid_xz.pdf");

  }
  //_____________ END TEMPORARY PLOTS


 //start a Timer
  //TTimer *timer = new TTimer(20);
  //  timer->SetCommand("Animate()");
  // timer->TurnOn();   
  // gPad->GetView()->RotateView(-90,180);
  //  gPad->Modified();
  // gPad->Update();

   //c2->SaveAs("myalgo_cm_xy.pdf");

   /* TTimer *timer = new TTimer(5);
   timer->SetCommand("Animate()");
   timer->TurnOn();*/
   //sleep(5);
   /*   TTimer *timer2 = new TTimer(5);
	timer2->SetCommand("*/
}
