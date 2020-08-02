{
  // gStyle->SetCanvasPreferGL(kTRUE);
  //    TCanvas *glc = new TCanvas("glc", "", 400, 0, 900, 900);
  //   gROOT->LoadMacro("$VMCWORKDIR/gconfig/rootlogon.C");
  //   rootlogon();
  
  gROOT->Reset();
  TFile f("evtcomplete_par.root", "READ");
  // if (!gGeoManager) {
  //   fairbasepar = (FairBaseParSet*) parfile.Get("FairBaseParSet");
  //   geoMan = gGeoManager;
  // }

  FairBaseParSet *parbase = 0;
  parbase = (FairBaseParSet*) f.Get("FairGeoParSet"); 
  TGeoManager *geoMan = 0;
  geoMan = gGeoManager;


  if(!geoMan) { 
    cout << "Could not find valid GeoManager. Abort now!" << endl; 
    exit(1); 
  } 
  
  parbase->GetName();
  PndGeoHandling geoH(1583944737,"evtcomplete_par.root");
  geoH.FillSensorMap();

  //FairBaseParSet *parbase = f.Get("FairBaseParSet");

  PndGeoSttPar *parameters = (PndGeoSttPar*) f.Get("PndGeoSttPar");
  PndSttMapCreator *mapper = new PndSttMapCreator(parameters);
  TClonesArray *fTubeArray = mapper->FillTubeArray();
  
  TCanvas *c = new TCanvas("c", "", 300, 300);

  TH2F *h2 = new TH2F("h2", "", 100, -42, 42, 100, -42, 42);
  h2->SetStats(kFALSE);
  h2->Draw();
  h2->GetXaxis()->SetTitle("x [cm]");
  h2->GetYaxis()->SetTitle("y [cm]");
  
  TArc *arc = NULL;
  TMarker *mrk = NULL;
  PndSttTube *tube = NULL;
  c->Clear();
  c->Draw();
  h2->Draw();
  
  int maxNumNei = -1;
  int minNumNei = 10000;
  int minID = 0;
  int maxID = 0;

  // draw all the tubes
  for(int itube = 1; itube < fTubeArray->GetEntriesFast(); itube++) {
    tube = (PndSttTube*) fTubeArray->At(itube);
    if(tube->IsParallel()) {
      arc = new TArc(tube->GetPosition().X(), tube->GetPosition().Y(), 0.5);
      //arc->SetFillColor(8);
      arc->SetLineColor(14);//17, 8
      arc->Draw("SAME");
    }
    else {
      TVector3 direction = tube->GetWireDirection();
      mrk = new TMarker(tube->GetPosition().X(), tube->GetPosition().Y(), 2);
      mrk->SetMarkerSize(0.4);
      if(direction.Phi() > 0){
	mrk->SetMarkerColor(2);
      }
      else {
	mrk->SetMarkerColor(4);
      }
      mrk->Draw("SAME");
    }
    // Find minimum and max number of neighbours
    TArrayI neigList = tube->GetNeighborings();
    int NumNeigh = neigList.GetSize();
    //Find max
    if( NumNeigh > maxNumNei ) {
      maxNumNei = NumNeigh;
      maxID = tube->GetTubeID();
    }
    if( NumNeigh < minNumNei ) {
      minNumNei = NumNeigh;
      minID = tube->GetTubeID();
    }
    
  }// END Tubes loop
  c->Update();
  c->Modified();
  // Save to pdf file
  c->SaveAs("STTGrid_XYPlane.pdf");

  // Print min and max num neighbours
  std::cout << "Min num neighbours = " << minNumNei << " tubeID = " << minID
	    << " Max num Neighb = " << maxNumNei << " tubeID = " << maxID
	    << '\n';
}
