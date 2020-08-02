/*************************************
 * Author: M. Babai (M.Babai@rug.nl) *
 * Version:                          *
 * License:                          *
 *************************************/
// #include <iostream>
// ///////////////////////
// #include "TNtuple.h"
// #include "TFile.h"
// #include "TCanvas.h"
// #include "TH2.h"

#define SAVE_PDF_PLOT 1

void PlotErrorValues(std::string const& Infile = "Tracks_Norm_1000Evts_08Gev_3D.root")
{
  Double_t const ww = 400;
  Double_t const hh = 300;
  
  TFile inp(Infile.c_str(),"READ");
  TNtuple *Error = (TNtuple*) inp.Get("ErrorEstimate");

  // Create plots
  // Undermerge Error
  TCanvas *c1 = new TCanvas("c1", "ErrorUnder", 400, 300 );
  Error->Draw("Error_underSegmentNorm >> hu(100,0,0.6");
  hu->SetTitle("Under-Merge");
  hu->SetLineColor(2);
  hu->Draw();
  // c1->SetLogy();
  //c1->SaveAs("undermerge08GevMu.eps");

  // Over Merge error
  TCanvas *c2 = new TCanvas("c2", "ErrorOver",  400, 300 );
  Error->Draw("Error_overSegmentNorm >> ho(100,0,0.6");
  ho->SetTitle("Over-Merge");
  ho->SetLineColor(2);
  // ho->Draw();
  // c2->SetLogy();
  // c2->Update()
  //c2->SaveAs("overmerge08GevMu.eps");

  // Total Error
  TCanvas *c3 = new TCanvas("c3", "ErrorTotal", 400, 300 );
  Error->Draw("TotalErrorNorm >> ht(100,0,0.6");
  ht->SetTitle("Total Error");
  ht->SetLineColor(2);
  // ht->Draw();
  // c3->SetLogy();
  //c3->SaveAs("totalError08GevMu.eps");
}
