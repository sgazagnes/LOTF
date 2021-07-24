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
  Double_t const ww = 500;
  Double_t const hh = 500;
  
  TFile inp(Infile1.c_str(),"READ");
  TNtuple *Error = (TNtuple*) inp.Get("ErrorEstimate");

  TNtuple *ErrorPT = (TNtuple*) inp.Get("PerTrackError");

  TNtuple *CurvPT = (TNtuple*) inp.Get("PerTrackCurv");

  TFile inp2(Infile2.c_str(),"READ");
  TNtuple *ErrorB = (TNtuple*) inp2.Get("ErrorEstimate");

  TNtuple *ErrorPTB = (TNtuple*) inp2.Get("PerTrackError");
  
  TNtuple *CurvPTB = (TNtuple*) inp2.Get("PerTrackCurv");
  // Create plots
  // Undermerge Error
  
  {
    TCanvas *c1 = new TCanvas("c1", "ErrorUnder", ww, hh );

    Error->SetLineColor(4);
    // Error->SetLineWidth(3);
    Error->Draw("Error_underMergeNorm");

    TH2F *hu =  (TH2F*)gPad->GetPrimitive("htemp");
    // hu->SetTitle("Under-Merge");
    hu->GetYaxis()->SetLimits(0,5);
    hu->GetXaxis()->SetLimits(-0.1,1.1);

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
    c1->SaveAs("Tot_Error_undermerge.png");

  }
  
  {
    // Over Merge error
    TCanvas *c2 = new TCanvas("c2", "ErrorOver",   ww, hh );

    Error->SetLineColor(4);
    Error->SetLineWidth(2);
    Error->Draw("Error_overMergeNorm");
    TH2F *ho = (TH2F*)gPad->GetPrimitive("htemp");

    ho->SetTitle("Over-Merge");
    ho->GetYaxis()->SetLimits(0,5);
    ho->GetXaxis()->SetLimits(-0.1,1.1);
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
    c2->SaveAs("Tot_Error_overmerge.png");

  }
  // Total Error
  {
    TCanvas *c3 = new TCanvas("c3", "ErrorTotal", ww, hh );
    Error->SetLineColor(4);
    Error->SetLineWidth(2);
    Error->Draw("TotalErrorNorm");
    TH2F *ht = (TH2F*)gPad->GetPrimitive("htemp");
    ht->SetTitle("Total Error");
    ht->GetYaxis()->SetLimits(0,5);
    ht->GetXaxis()->SetLimits(-0.1,1.1);
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
    c3->SaveAs("Tot_Error.png");

    //c3->SaveA    c1->SaveAs("Tot_Error_undermerge.png");
;
  }
  {
    TCanvas *c4 = new TCanvas("c4", "JacardSingle",  ww, hh);
    ErrorPT->SetLineColor(4);
    ErrorPT->SetLineWidth(2);
    ErrorPT->Draw("Jacardsingle");
    TH2F *hp = (TH2F*)gPad->GetPrimitive("htemp");

    hp->SetTitle("Jacardsingle");
    hp->GetYaxis()->SetLimits(0,5);
    hp->GetXaxis()->SetLimits(-0.1,1.1);

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
 
    //c4->Update();
    c4->SaveAs("Jacardsin.png");
    }
  // Over Merge error
   {
    TCanvas *c5 = new TCanvas("c5", "F1 score",   ww, hh );
    ErrorPT->SetLineColor(4);
    // ErrorPT->SetLineWidth(3);
    ErrorPT->Draw("Jacardaverage");
    TH2F *ho = (TH2F*)gPad->GetPrimitive("htemp");

    ho->SetTitle("f1score");
    ho->GetYaxis()->SetLimits(0,5);
    ho->GetXaxis()->SetLimits(-0.1,1.1);
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
    c5->SaveAs("f1score.png");

    }
  /*{// Total Error
    TCanvas *c6 = new TCanvas("c6", "Found",  ww, hh );
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
    c6->SaveAs("Found.png");
    }*/
     {
    TCanvas *c7 = new TCanvas("c7", "deltaRCanvas",   ww, hh );
    CurvPT->SetLineColor(4);
    // ErrorPT->SetLineWidth(3);
    CurvPT->Draw("(((1/MC_r)-(1/tr_r))/(1/MC_r))");
    TH2F *ho = (TH2F*)gPad->GetPrimitive("htemp");

    ho->SetTitle("deltaRCanvas");
    ho->SetXTitle("(R_{MC} - R_{Tr}) [cm]");
    ho->SetYTitle("counts");
    //   ho->GetYaxis()->SetLimits(0,5);
    ho->GetXaxis()->SetLimits(-50,50);
    c7->Modified();

     c7->Update();
     TPaveStats *stats =(TPaveStats*)c7->GetPrimitive("stats");
    stats->SetName("h1stats");
    stats->SetY1NDC(.9);
    stats->SetY2NDC(.7);
    stats->SetTextColor(4);
     CurvPTB->SetLineColor(2);
    CurvPTB->Draw("(((1/MC_r)-(1/tr_r))/(1/MC_r))", "","sames");
    c7->Update();
    TPaveStats *stats2 = (TPaveStats*)c7->GetPrimitive("stats");
    stats2->SetName("h1stats2");
    stats2->SetY1NDC(.65);
    stats2->SetY2NDC(.45);
    stats2->SetTextColor(2);
    c7->Update();
    //   c7->SaveAs("Jacarave.png");


    }
  {// Total Error
    TCanvas *c9 = new TCanvas("c9", "Displacement X",ww, hh );
    ErrorPT->SetLineColor(4);
    // ErrorPT->SetLineWidth(3);
    ErrorPT->Draw("disX");
    TH2F *ht = (TH2F*)gPad->GetPrimitive("htemp");

    ht->SetTitle("Displacement X");
    ht->GetYaxis()->SetLimits(0,5);
    ht->GetXaxis()->SetLimits(-0.1,2);
    c9->Modified();

    c9->Update();
    TPaveStats *stats =(TPaveStats*)c9->GetPrimitive("stats");
    stats->SetName("h1stats");
    stats->SetY1NDC(.9);
    stats->SetY2NDC(.7);
    stats->SetTextColor(4);
    ErrorPTB->SetLineColor(2);
    ErrorPTB->Draw("disX", "","sames");
    c9->Update();
    TPaveStats *stats2 = (TPaveStats*)c9->GetPrimitive("stats");
    stats2->SetName("h1stats2");
    stats2->SetY1NDC(.65);
    stats2->SetY2NDC(.45);
    stats2->SetTextColor(2);
    c9->Update();
    c9->SaveAs("Disx.png");

  }
        {// T10tal Error
    TCanvas *c10 = new TCanvas("c10", "Displacement Y",  ww, hh );
    ErrorPT->SetLineColor(4);
    // ErrorPT->SetLineWidth(3);
    ErrorPT->Draw("disY");
    TH2F *ht = (TH2F*)gPad->GetPrimitive("htemp");

    ht->SetTitle("Displacement Y");
    ht->GetYaxis()->SetLimits(0,5);
    ht->GetXaxis()->SetLimits(-0.1,2);
    c10->Modified();

    c10->Update();
    TPaveStats *stats =(TPaveStats*)c10->GetPrimitive("stats");
    stats->SetName("h1stats");
    stats->SetY1NDC(.9);
    stats->SetY2NDC(.7);
    stats->SetTextColor(4);
    ErrorPTB->SetLineColor(2);
    ErrorPTB->Draw("disY", "","sames");
    c10->Update();
    TPaveStats *stats2 = (TPaveStats*)c10->GetPrimitive("stats");
    stats2->SetName("h1stats2");
    stats2->SetY1NDC(.65);
    stats2->SetY2NDC(.45);
    stats2->SetTextColor(2);
    c10->Update();
    c10->SaveAs("Disy.png");

  }
	    {// Total Error
    TCanvas *c11 = new TCanvas("c11", "Displacement Z",  ww, hh);
    ErrorPT->SetLineColor(4);
    // ErrorPT->SetLineWidth(3);
    ErrorPT->Draw("disZ");
    TH2F *ht = (TH2F*)gPad->GetPrimitive("htemp");

    ht->SetTitle("Displacement Z");
        ht->GetYaxis()->SetLimits(0,5);
    ht->GetXaxis()->SetLimits(-0.1,40);
    c11->Modified();

    c11->Update();
    TPaveStats *stats =(TPaveStats*)c11->GetPrimitive("stats");
    stats->SetName("h1stats");
    stats->SetY1NDC(.9);
    stats->SetY2NDC(.7);
    stats->SetTextColor(4);
    ErrorPTB->SetLineColor(2);
    ErrorPTB->Draw("disZ", "","sames");
    c11->Update();
    TPaveStats *stats2 = (TPaveStats*)c11->GetPrimitive("stats");
    stats2->SetName("h1stats2");
    stats2->SetY1NDC(.65);
    stats2->SetY2NDC(.45);
    stats2->SetTextColor(2);
    c11->Update();
    c11->SaveAs("Disz.png");

    }

  





  /*COMPLEX*//*
  {
    TCanvas *c1 = new TCanvas("c1", "ErrorUnder", ww, hh );
 
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

    c1->Update();
    c1->SaveAs("Tot_Error_undermerge.png");

  }
  
  {
    // Over Merge error
    TCanvas *c2 = new TCanvas("c2", "ErrorOver",   ww, hh );

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
    c2->SaveAs("Tot_Error_overmerge.png");

  }
  // Total Error
  {
    TCanvas *c3 = new TCanvas("c3", "ErrorTotal", ww, hh );
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
    c3->SaveAs("Tot_Error.png");

    //c3->SaveA    c1->SaveAs("Tot_Error_undermerge.png");
;
  }
  {
    TCanvas *c4 = new TCanvas("c4", "JacardSingle",  ww, hh);
    ErrorPT->SetLineColor(4);
    ErrorPT->SetLineWidth(2);
    ErrorPT->Draw("Jacardsingle","(complex == 1)");
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
    ErrorPTB->Draw("Jacardsingle", "(complex == 1)","sames");
    c4->Update();
    TPaveStats *stats2 = (TPaveStats*)c4->GetPrimitive("stats");
    stats2->SetName("h1stats2");
    stats2->SetY1NDC(.65);
    stats2->SetY2NDC(.45);
    stats2->SetTextColor(2);
    c4->Update();
    c4->SaveAs("Jacardsin.png");
   }
  // Over Merge error
  /* {
    TCanvas *c5 = new TCanvas("c5", "JacardAverage",   ww, hh );
    ErrorPT->SetLineColor(4);
    // ErrorPT->SetLineWidth(3);
    ErrorPT->Draw("Jacardaverage","(complex == 1)");
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
    ErrorPTB->Draw("Jacardaverage", "(complex == 1)","sames");
    c5->Update();
    TPaveStats *stats2 = (TPaveStats*)c5->GetPrimitive("stats");
    stats2->SetName("h1stats2");
    stats2->SetY1NDC(.65);
    stats2->SetY2NDC(.45);
    stats2->SetTextColor(2);
    c5->Update();
    c5->SaveAs("Jacarave.png");

    }
  {// Total Error
    TCanvas *c6 = new TCanvas("c6", "Found",  ww, hh );
    ErrorPT->SetLineColor(4);
    //ErrorPT->SetLineWidth(3);
    ErrorPT->Draw("misMatched", "(complex == 1)");
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
    ErrorPTB->Draw("misMatched", "(complex == 1)","sames");
    c6->Update();
    TPaveStats *stats2 = (TPaveStats*)c6->GetPrimitive("stats");
    stats2->SetName("h1stats2");
    stats2->SetY1NDC(.65);
    stats2->SetY2NDC(.45);
    stats2->SetTextColor(2);
    c6->Update();
    c6->SaveAs("Found.png");
  }
  {
    TCanvas *c7 = new TCanvas("c7", "UnderMergeperTrack",  ww, hh);
    ErrorPT->SetLineColor(4);
    // ErrorPT->SetLineWidth(3);
    ErrorPT->Draw("UnderMergeError", "(complex == 1)");
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
    ErrorPTB->Draw("UnderMergeError", "(complex == 1)","sames");
    c7->Update();
    TPaveStats *stats2 = (TPaveStats*)c7->GetPrimitive("stats");
    stats2->SetName("h1stats2");
    stats2->SetY1NDC(.65);
    stats2->SetY2NDC(.45);
    stats2->SetTextColor(2);
    c7->Update();
    c7->SaveAs("Undermerge_track.png");

  }
  {// Total Error
    TCanvas *c8 = new TCanvas("c8", "Overmegrepertrack",  ww, hh );
    ErrorPT->SetLineColor(4);
    // ErrorPT->SetLineWidth(3);
    ErrorPT->Draw("OverMergeError", "(complex == 1)");
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
    ErrorPTB->Draw("OverMergeError", "(complex == 1)","sames");
    c8->Update();
    TPaveStats *stats2 = (TPaveStats*)c8->GetPrimitive("stats");
    stats2->SetName("h1stats2");
    stats2->SetY1NDC(.65);
    stats2->SetY2NDC(.45);
    stats2->SetTextColor(2);
    c8->Update();
    c8->SaveAs("Overmerge_track.png");

  }
    {// Total Error
    TCanvas *c9 = new TCanvas("c9", "Displacement X",ww, hh );
    ErrorPT->SetLineColor(4);
    // ErrorPT->SetLineWidth(3);
    ErrorPT->Draw("disX", "(complex == 1)");
    TH2F *ht = (TH2F*)gPad->GetPrimitive("htemp");

    ht->SetTitle("Displacement X");
    c9->Modified();

    c9->Update();
    TPaveStats *stats =(TPaveStats*)c9->GetPrimitive("stats");
    stats->SetName("h1stats");
    stats->SetY1NDC(.9);
    stats->SetY2NDC(.7);
    stats->SetTextColor(4);
    ErrorPTB->SetLineColor(2);
    ErrorPTB->Draw("disX", "(complex == 1)","sames");
    c9->Update();
    TPaveStats *stats2 = (TPaveStats*)c9->GetPrimitive("stats");
    stats2->SetName("h1stats2");
    stats2->SetY1NDC(.65);
    stats2->SetY2NDC(.45);
    stats2->SetTextColor(2);
    c9->Update();
    c9->SaveAs("Disx.png");

  }
        {// T10tal Error
    TCanvas *c10 = new TCanvas("c10", "Displacement Y",  ww, hh );
    ErrorPT->SetLineColor(4);
    // ErrorPT->SetLineWidth(3);
    ErrorPT->Draw("disY", "(complex == 1)");
    TH2F *ht = (TH2F*)gPad->GetPrimitive("htemp");

    ht->SetTitle("Displacement Y");
    c10->Modified();

    c10->Update();
    TPaveStats *stats =(TPaveStats*)c10->GetPrimitive("stats");
    stats->SetName("h1stats");
    stats->SetY1NDC(.9);
    stats->SetY2NDC(.7);
    stats->SetTextColor(4);
    ErrorPTB->SetLineColor(2);
    ErrorPTB->Draw("disY", "(complex == 1)","sames");
    c10->Update();
    TPaveStats *stats2 = (TPaveStats*)c10->GetPrimitive("stats");
    stats2->SetName("h1stats2");
    stats2->SetY1NDC(.65);
    stats2->SetY2NDC(.45);
    stats2->SetTextColor(2);
    c10->Update();
    c10->SaveAs("Disy.png");

  }
	    {// Total Error
    TCanvas *c11 = new TCanvas("c11", "Displacement Z",  ww, hh);
    ErrorPT->SetLineColor(4);
    // ErrorPT->SetLineWidth(3);
    ErrorPT->Draw("disZ", "(complex == 1)");
    TH2F *ht = (TH2F*)gPad->GetPrimitive("htemp");

    ht->SetTitle("Displacement Z");
    c11->Modified();

    c11->Update();
    TPaveStats *stats =(TPaveStats*)c11->GetPrimitive("stats");
    stats->SetName("h1stats");
    stats->SetY1NDC(.9);
    stats->SetY2NDC(.7);
    stats->SetTextColor(4);
    ErrorPTB->SetLineColor(2);
    ErrorPTB->Draw("disZ", "(complex == 1)","sames");
    c11->Update();
    TPaveStats *stats2 = (TPaveStats*)c11->GetPrimitive("stats");
    stats2->SetName("h1stats2");
    stats2->SetY1NDC(.65);
    stats2->SetY2NDC(.45);
    stats2->SetTextColor(2);
    c11->Update();
    c11->SaveAs("Disz.png");

    }*/
}
