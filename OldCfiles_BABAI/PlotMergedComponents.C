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

#define SAVE_PDF_PLOT 1

void PlotMergedComponents( std::string const& Infile = "Track3DPositions.root",
			   Double_t const ww    = 990,
			   Double_t const hh    = 990,
			   size_t   const nEv   = 1,
			   size_t   const nComp = 9,
			   size_t   const dim   = 3,
			   bool     const showGrid = false
			   )
{
  TFile inp(Infile.c_str(),"READ");
  std::string tpName;
  std::vector< std::vector<std::string>* >AllComponents;
  
  // Build a list of aal ntuples containing the data for different
  // orientations.
  // Event loop
  for(size_t ev = 0; ev < nEv; ++ev) {
    std::vector<std::string>* tNames = new std::vector<std::string>();
    std::stringstream evN;
    evN << ev;
    //Orientations
    for(size_t i = 0; i < nComp; i++) {
      std::stringstream oriN;
      oriN << i;
      tpName = "MergedComp_Evt_" + evN.str() + "_" + oriN.str();
      //std::cout << tpName << std::endl;
      tNames->push_back(tpName);
    }//End Of orient. loop
    AllComponents.push_back(tNames);
  }//END event loop

   // Create NTuples and plot
  TCanvas *compCanvas = new TCanvas("AllMergedComponentsCanvas", "AllMergedComponents", ww, hh );
  
  TNtuple* grid = (TNtuple*) inp.Get("ExtendedGrid");
  if(!showGrid) {
    grid->SetMarkerColor(kWhite);
  }
  if(dim == 2) {
    grid->Draw("y:x","","");
  }
  else {
    grid->Draw("y:x:z","","");
  }

  for(size_t i = 0; i < AllComponents.size(); ++i) {
    std::stringstream pltCnt;
    pltCnt << i;
    std::string fname = "plot_" + pltCnt.str() + ".pdf";
    // Fetch the list of plots for the current event.
    std::vector<std::string> const *currEvt = AllComponents[i];
    int color = 1;
    for(size_t j = 0; j < currEvt->size(); ++j) {
      std::string const &name = currEvt->at(j);
      TNtuple* plot = (TNtuple*) inp.Get(name.c_str());
      // color = color + (j%9);
      color += j;
      // if(color == 5) {
      if( (j % 4) == 0) {
	color = kBlack;
      }
      if( (j % 4) == 1) {
	color =  kRed ;
      }
      if( (j % 4) == 2) {
	color = kBlue;
      }
      if( (j % 4) == 3) {
	color = kGreen + 2;
      }
      plot->SetMarkerColor(color);
      plot->SetMarkerStyle( (j%30) + 20);
      plot->SetMarkerSize(0.8);
      if(dim == 2) {
	plot->Draw("y:x","","same");
      }
      if(dim == 3){
	plot->Draw("y:x:z_Det","","same");
      }
    }
#if ( SAVE_PDF_PLOT > 1 )
    c1->SaveAs(fname.c_str());
#endif
    compCanvas->Update();
  }
  //exit(0);
}
