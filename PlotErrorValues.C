/*************************************
 * Author: M. Babai (M.Babai@rug.nl) *
 * Version:                          *
 * License:                          *
 *************************************/
// #include <iostream>
// ///////////////////////

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

#define SAVE_PDF_PLOT 1

void PlotErrorValues(std::string const& Infile1 = "Tracks_output.root",std::string const& Infile2 = "../Test_babai/Tracks_output.root")
{
  Double_t const ww = 400;
  Double_t const hh = 400;
  
  TFile inp(Infile1.c_str(),"READ");
  TNtuple *Error = (TNtuple*) inp.Get("ErrorEstimate");

  TNtuple *ErrorPT = (TNtuple*) inp.Get("PerTrackError");

  TFile inp2(Infile2.c_str(),"READ");
  TNtuple *ErrorB = (TNtuple*) inp2.Get("ErrorEstimate");

  TNtuple *ErrorPTB = (TNtuple*) inp2.Get("PerTrackError");
  // Create plots
  // Undermerge Error
  {
    TCanvas *c1 = new TCanvas("c1", "ErrorUnder", 400, 300 );
 
    Error->SetLineColor(4);
    // Error->SetLineWidth(3);
    Error->Draw("Error_underMergeNorm");
    //  TH2F *hu =  (TH2F*)gPad->GetPrimitive("htemp");
    // hu->SetTitle("Under-Merge");


    // TLegend  *legend = new TLegend(0.80,0.3,0.995,0.4); // we need different positions for the legend to not 

    //legend->AddEntry(Error,"ONE","L");
    c1->Modified();

    c1->Update();
    TPaveStats *stats =(TPaveStats*)c1->GetPrimitive("stats");
    stats->SetName("h1stats");
    stats->SetY1NDC(.9);
    stats->SetY2NDC(.7);
    stats->SetTextColor(4);
    //  TList *list = stats->GetListOfLines();
    // TText *title_line = new TText(0,0,"Stats Box Title");
    //  list->AddFirst(title_line); //add new line to top of stats box ??
    //TH2F *hu2 = (TH2F*)gPad->GetPrimitive("h2");

    //TPaveStats stats1 = hu->GetListOfFunctions().FindObject("stats").Clone("stats1");
    //stats1.SetY1NDC(.5);
    //stats1.SetY2NDC(.7);
    //hu->SetLineColor(2);
    ErrorB->SetLineColor(2);
    ErrorB->Draw("Error_underMergeNorm", "","sames");
    // TH2F *hb =  (TH2F*)gPad->GetPrimitive("h2");
    //hu->SetTitle("Under-Merge");
    c1->Update();
    TPaveStats *stats2 = (TPaveStats*)c1->GetPrimitive("stats");
    stats2->SetName("h1stats2");
    stats2->SetY1NDC(.65);
    stats2->SetY2NDC(.45);
    stats2->SetTextColor(2);
    // legend->AddEntry(ErrorB,"TWO","L");
    // legend->Draw();
    // c1->SetLogy();
    //c1->SaveAs("undermerge08GevMu.eps");
    stats->Draw();
    stats2->Draw();
    c1->Update();
  }
  
  {
    // Over Merge error
    TCanvas *c2 = new TCanvas("c2", "ErrorOver",  400, 300 );

    Error->SetLineColor(4);
    Error->SetLineWidth(2);
    Error->Draw("Error_overMergeNorm");
    TH2F *ho = (TH2F*)gPad->GetPrimitive("htemp");

    ho->SetTitle("Over-Merge");

    c2->Modified();

    c2->Update();
    TPaveStats *stats =(TPaveStats*)c2->GetPrimitive("stats");
    stats->SetName("h1stats");
    stats->SetY1NDC(.9);
    stats->SetY2NDC(.7);
    stats->SetTextColor(4);
  
    ErrorB->SetLineColor(2);
    ErrorB->Draw("Error_overMergeNorm", "","sames");
    c2->Update();
    TPaveStats *stats2 = (TPaveStats*)c2->GetPrimitive("stats");
    stats2->SetName("h1stats2");
    stats2->SetY1NDC(.65);
    stats2->SetY2NDC(.45);
    stats2->SetTextColor(2);
    // c1->SetLogy();
    //c1->SaveAs("undermerge08GevMu.eps");
    c2->Update();
  }
  // Total Error
  {
    TCanvas *c3 = new TCanvas("c3", "ErrorTotal", 400, 300 );
    Error->SetLineColor(4);
    Error->SetLineWidth(2);
    Error->Draw("TotalErrorNorm");
    TH2F *ht = (TH2F*)gPad->GetPrimitive("htemp");
    ht->SetTitle("Total Error");

    c3->Modified();

    c3->Update();
    TPaveStats *stats =(TPaveStats*)c3->GetPrimitive("stats");
    stats->SetName("h1stats");
    stats->SetY1NDC(.9);
    stats->SetY2NDC(.7);
    stats->SetTextColor(4);
    ErrorB->SetLineColor(2);
    ErrorB->Draw("TotalErrorNorm", "","sames");
    c3->Update();
    TPaveStats *stats2 = (TPaveStats*)c3->GetPrimitive("stats");
    stats2->SetName("h1stats2");
    stats2->SetY1NDC(.65);
    stats2->SetY2NDC(.45);
    stats2->SetTextColor(2);
    c3->Update();

    //c3->SaveAs("totalError08GevMu.eps");
  }
  {
    TCanvas *c4 = new TCanvas("c4", "JacardSingle", 400, 300 );
    ErrorPT->SetLineColor(4);
    ErrorPT->SetLineWidth(2);
    ErrorPT->Draw("Jacardsingle");
    TH2F *hu = (TH2F*)gPad->GetPrimitive("htemp");


    hu->SetTitle("Jacardsingle");
    c4->Modified();

    c4->Update();
    TPaveStats *stats =(TPaveStats*)c4->GetPrimitive("stats");
    stats->SetName("h1stats");
    stats->SetY1NDC(.9);
    stats->SetY2NDC(.7);
    stats->SetTextColor(4);
    ErrorPTB->SetLineColor(2);
    ErrorPTB->Draw("Jacardsingle", "","sames");
    c4->Update();
    TPaveStats *stats2 = (TPaveStats*)c4->GetPrimitive("stats");
    stats2->SetName("h1stats2");
    stats2->SetY1NDC(.65);
    stats2->SetY2NDC(.45);
    stats2->SetTextColor(2);
    c4->Update();

  }
  // Over Merge error
  {
    TCanvas *c5 = new TCanvas("c5", "JacardAverage",  400, 300 );
    ErrorPT->SetLineColor(4);
    // ErrorPT->SetLineWidth(3);
    ErrorPT->Draw("Jacardaverage");
    TH2F *ho = (TH2F*)gPad->GetPrimitive("htemp");

    ho->SetTitle("Jacardaverae");
    c5->Modified();

    c5->Update();
    TPaveStats *stats =(TPaveStats*)c5->GetPrimitive("stats");
    stats->SetName("h1stats");
    stats->SetY1NDC(.9);
    stats->SetY2NDC(.7);
    stats->SetTextColor(4);
    ErrorPTB->SetLineColor(2);
    ErrorPTB->Draw("Jacardaverage", "","sames");
    c5->Update();
    TPaveStats *stats2 = (TPaveStats*)c5->GetPrimitive("stats");
    stats2->SetName("h1stats2");
    stats2->SetY1NDC(.65);
    stats2->SetY2NDC(.45);
    stats2->SetTextColor(2);
    c5->Update();

  }
  {// Total Error
    TCanvas *c6 = new TCanvas("c6", "Found", 400, 300 );
    ErrorPT->SetLineColor(4);
    //ErrorPT->SetLineWidth(3);
    ErrorPT->Draw("misMatched");
    TH2F *ht = (TH2F*)gPad->GetPrimitive("htemp");

    ht->SetTitle("Total Error");
    c6->Modified();

    c6->Update();
    TPaveStats *stats =(TPaveStats*)c6->GetPrimitive("stats");
    stats->SetName("h1stats");
    stats->SetY1NDC(.9);
    stats->SetY2NDC(.7);
    stats->SetTextColor(4);
    ErrorPTB->SetLineColor(2);
    ErrorPTB->Draw("misMatched", "","sames");
    c6->Update();
    TPaveStats *stats2 = (TPaveStats*)c6->GetPrimitive("stats");
    stats2->SetName("h1stats2");
    stats2->SetY1NDC(.65);
    stats2->SetY2NDC(.45);
    stats2->SetTextColor(2);
    c6->Update();

  }
  {
    TCanvas *c7 = new TCanvas("c7", "UnderMergeperTrack",  400, 300 );
    ErrorPT->SetLineColor(4);
    // ErrorPT->SetLineWidth(3);
    ErrorPT->Draw("UnderMergeError");
    TH2F *ho = (TH2F*)gPad->GetPrimitive("htemp");

    ho->SetTitle("UnderMergeError");
    c7->Modified();

    c7->Update();
    TPaveStats *stats =(TPaveStats*)c7->GetPrimitive("stats");
    stats->SetName("h1stats");
    stats->SetY1NDC(.9);
    stats->SetY2NDC(.7);
    stats->SetTextColor(4);
    ErrorPTB->SetLineColor(2);
    ErrorPTB->Draw("UnderMergeError", "","sames");
    c7->Update();
    TPaveStats *stats2 = (TPaveStats*)c7->GetPrimitive("stats");
    stats2->SetName("h1stats2");
    stats2->SetY1NDC(.65);
    stats2->SetY2NDC(.45);
    stats2->SetTextColor(2);
    c7->Update();

  }
  {// Total Error
    TCanvas *c8 = new TCanvas("c8", "Overmegrepertrack", 400, 300 );
    ErrorPT->SetLineColor(4);
    // ErrorPT->SetLineWidth(3);
    ErrorPT->Draw("OverMergeError");
    TH2F *ht = (TH2F*)gPad->GetPrimitive("htemp");

    ht->SetTitle("OverMergeError");
    c8->Modified();

    c8->Update();
    TPaveStats *stats =(TPaveStats*)c8->GetPrimitive("stats");
    stats->SetName("h1stats");
    stats->SetY1NDC(.9);
    stats->SetY2NDC(.7);
    stats->SetTextColor(4);
    ErrorPTB->SetLineColor(2);
    ErrorPTB->Draw("OverMergeError", "","sames");
    c8->Update();
    TPaveStats *stats2 = (TPaveStats*)c8->GetPrimitive("stats");
    stats2->SetName("h1stats2");
    stats2->SetY1NDC(.65);
    stats2->SetY2NDC(.45);
    stats2->SetTextColor(2);
    c8->Update();

  }
}
