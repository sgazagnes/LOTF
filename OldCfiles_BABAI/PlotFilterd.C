/*************************************
 * Author: M. Babai (M.Babai@rug.nl) *
 * Version:                          *
 * License:                          *
 *************************************/
#include <string>
#include <sstream>
#include <vector>


#include "TNtuple.h"
#include "TFile.h"

#define WRITE_PLOTS 0
#define SAVE_PDF_PLOT 0

void PlotFilterd( std::string const& InfileAll = "Track3DPositions.root",
                  Double_t const ww = 600,
                  Double_t const hh = 600
                  )
{
  TFile inpAll(InfileAll.c_str(), "READ");

  // Fetch data
  TNtuple* PosOrig    = (TNtuple*) inpAll.Get("Pos");
  TNtuple* PosFilterd = (TNtuple*) inpAll.Get("Coord");

  // _________________ Test multiple ______________
  TCanvas *c1 = new TCanvas("Filterd", "Coordinates All", (ww + 600), (hh + 00) );
  TLatex Tl;
  Tl.SetTextSize(0.03);

  c1->Divide(2,1);
  // Original coordinates.
  c1->cd(1);
  PosOrig->SetMarkerColor(4);
  PosOrig->SetMarkerStyle(6);
  PosOrig->Draw("y:x","","");
  Tl.DrawLatex(-10.0,45.0,"Original");

  // Filterd coordinates.
  c1->cd(2);
  PosFilterd->SetMarkerColor(2);
  PosFilterd->SetMarkerStyle(6);
  PosFilterd->Draw("y:x","","");
  //Tl.DrawLatex(-10.0,45.0,"Area Opening #lambda = 12");

  //c1->SaveAs("Area_Open_Lambda12.eps");
}
