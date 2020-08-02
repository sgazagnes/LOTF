#include "TFile.h"
#include "TCanvas.h"
#include "TNtuple.h"

void DrawCorrected(char const *InFile = "Tracks_output.root")
{
  TCanvas *canvas = new TCanvas("CorrectedSkewed","CorrectedSkewed", 450, 450);
  
  TFile f(InFile,"READ");

  TNtuple *extendedGrid = (TNtuple*) f.Get("ExtendedGrid");

  extendedGrid->Draw("y:x","","");
  
  TNtuple *CollectedCoords = (TNtuple*) f.Get("CoordCollected");

  CollectedCoords->SetMarkerStyle(6);
  CollectedCoords->SetMarkerColor(4);
  
  CollectedCoords->Draw("y_Det:x_Det","", "same");

  canvas->Update();
}
