/*************************************
 * Author: M. Babai (M.Babai@rug.nl) *
 * Version:                          *
 * License:                          *
 *************************************/
#include <string>
#include <iostream>

#include "TNtuple.h"
#include "TCanvas.h"
#include "TFile.h"

#define WRITE_PLOTS 0
#define SAVE_PDF_PLOT 0

void PlotTracks( std::string const& Infile = "Tracks_output.root",
                 Double_t const ww = 500,
                 Double_t const hh = 500,
                 size_t   const nEv = 5,
		 size_t   const dim = 2
                 )
{
  std::cout << "Dummy Print PlotTracks" << nEv << '\n';
  TFile inp(Infile.c_str(),"READ");

  TNtuple* grid = (TNtuple*) inp.Get("ExtendedGrid");
  
  // Cartesian coordinates
  TNtuple* Pos   = (TNtuple*) inp.Get("Pos");
  TNtuple* MCpos = (TNtuple*) inp.Get("MCpos");

  // Polar coordinates.
  TNtuple* MCposPolar = (TNtuple*) inp.Get("MCposPolar");
  TNtuple* PosPolar   = (TNtuple*) inp.Get("PosPolar");

  Pos->AddFriend("MCpos");

  PosPolar->AddFriend("MCposPolar");

  // _________________ Test multiple ______________
  //for(size_t i = 0; i < nEv; ++i) {
  TCanvas *c1 = new TCanvas("EventPlot", "All Collected Coordinates", ww, hh );
  c1->Divide(2,2);

  //___ Cartesian coordinates
  c1->cd(1);
  Pos->SetMarkerColor(4);
  Pos->SetMarkerStyle(6);
  // Pos->SetMarkerSize(0.7);
  if(dim == 2) {
    grid->Draw("y:x","","");
    Pos->Draw("y:x","mcy > -10000","same");
  }
  if(dim == 3) {
    grid->Draw("y:x:z","","");
    Pos->Draw("y:x:z","mcy > -10000","same");
  }
  // MC data
  c1->cd(2);
  Pos->SetMarkerColor(2);
  //Pos->SetMarkerStyle(4);
  Pos->SetMarkerStyle(6);
  //Pos->SetMarkerSize(0.6);
  if(dim == 2) {
    grid->Draw("y:x","","");
    Pos->Draw("mcy:mcx","mcy > -10000","same");
  }
  if(dim == 3) {
    grid->Draw("y:x:z","","");
    Pos->Draw("mcy:mcx:mcz","mcy > -10000","same");
  }
  //_________________ Polar coordinates
  c1->cd(3);
  PosPolar->SetMarkerColor(4);
  PosPolar->SetMarkerStyle(6);
  //PosPolar->SetMarkerSize(0.7);
  PosPolar->Draw("thetaDeg:r","mcr > -10000","");

  // MC data
  c1->cd(4);
  PosPolar->SetMarkerColor(2);
  //PosPolar->SetMarkerStyle(4);
  PosPolar->SetMarkerStyle(6);
  //PosPolar->SetMarkerSize(0.6);
  PosPolar->Draw("mcthetaDeg:mcr","mcr > -10000","");
  // _________________ Test multiple ______________

#if ( WRITE_PLOTS > 0 )
  TFile pl("plots_tracks.root","RECREATE");
  c1->Write();
  pl.Close();
#endif

#if ( SAVE_PDF_PLOT > 1 )
  //c1->SaveAs("CollectionAllCoords.pdf");
#endif

#if ( SAVE_PDF_PLOT > 0 )
  //c1->SaveAs("CollectionAllCoords.pdf");
#endif
}
