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

  PndGeoSttPar *parameters =  (PndGeoSttPar*) f.Get("PndGeoSttPar");// 

  PndSttMapCreator *mapper = new PndSttMapCreator(parameters);

  TClonesArray *fTubeArray = mapper->FillTubeArray();


  TCanvas *c = new TCanvas("c", "", 400, 0, 900, 900);
  TH2F *h2 = new TH2F("h2", "", 100, -42, 42, 100, -42, 42);
  h2->SetStats(kFALSE);
  h2->Draw();
  
  TArc *arc = NULL;
  TMarker *mrk = NULL;
  PndSttTube *tube = NULL;
  PndSttTube *sttTube = nullptr;
  int nodeID;
  int nodeType;
  int SectorLimit = 0;
  int Sector = 0;
  int Layer = 0;
  bool LayerLimit = 0 ;
  int itube = 0;

  while(itube != -1) {
    c->Clear();
    c->Draw();
    h2->Draw();
    
    
    // draw all the tubes
    for(itube = 1; itube < fTubeArray->GetEntriesFast(); itube++) {
      tube = (PndSttTube*) fTubeArray->At(itube);
      //sttTube = PndSttTubeMap::Instance()->GetTube(iTube);

      if(tube->IsParallel()) {
	arc = new TArc(tube->GetPosition().X(), tube->GetPosition().Y(), 0.5);
	arc->Draw("SAME");
      }
      else {
	mrk = new TMarker(tube->GetPosition().X(), tube->GetPosition().Y(), 6);
	mrk->Draw("SAME");
      }

      Sector     = tube->GetSectorID();
      Layer      = tube->GetLayerID();
      LayerLimit = tube->IsLayerLimit();
      SectorLimit = tube->IsSectorLimit();
    }
    // ............................
    
    cout << "WHICH tube do you want to test? [to exit type -1]" << endl;
    cin >> itube;
    if(itube == -1) {
      exit(0);
    }
    cout << "TESTING tube " << itube << endl;
    tube = (PndSttTube*) fTubeArray->At(itube);
    if(!tube) {
      cout << "this tube does not exist" << endl; 
      //continue; 
      exit(1);
    }
    TVector3 pos = tube->GetPosition();
    Double_t hlfLength = tube->GetHalfLength();
    TArrayI neighborings = tube->GetNeighborings();
    std::cout << "Position X = " << pos.X() << " Y = " << pos.Y() << " Z = " << pos.Z()
	      << " Half Length = " << hlfLength
	      << " N neighbors: " << neighborings.GetSize()
	      << " Sector = " << tube->GetSectorID()
	      << " Layer = " << tube->GetLayerID()
	      << " Skewed = " << tube->IsSkew()
	      << " Parallel = " << tube->IsParallel()
	      << " Layer Limit = " << tube->IsLayerLimit()
	      << " Sector Limit = " <<  tube->IsSectorLimit()
	      << std::endl;
    TMatrixT<Double_t> mat = tube->GetRotationMatrix();
    mat.Print();
    
    std::cout << "Direction=\n";
    TVector3 direct = tube->GetWireDirection();
    direct.Print();
    
    arc = new TArc(tube->GetPosition().X(), tube->GetPosition().Y(), 0.5);
    arc->SetFillColor(kYellow);
    arc->Draw("SAME");
    
    for(int inei = 0; inei < neighborings.GetSize(); inei++) { 
      tube = (PndSttTube*) fTubeArray->At(neighborings.At(inei));
      TVector3 position = tube->GetPosition(); 
      cout << " " << neighborings.At(inei);
      if(tube->IsParallel()) {
	arc = new TArc(position.X(), position.Y(), 0.5); 
	arc->SetFillColor(3);
	arc->Draw("SAME");
      }
      else {
	//mrk = new TMarker(position.X(), position.Y(), 6);
	TVector3 extr1 = position - tube->GetWireDirection() * tube->GetHalfLength();
	TVector3 extr2 = position + tube->GetWireDirection() * tube->GetHalfLength();
        mrk = new TMarker(extr1.X(), extr1.Y(), 6);
        //mrk = new TMarker(extr2.X(), extr2.Y(), 6);
	mrk->SetMarkerColor(3);
	mrk->Draw("SAME");
      }
    }
    cout << endl;
    c->Update();
    c->Modified();
    //c->SaveAs("Tube_1053_Image.pdf");
  }
}
