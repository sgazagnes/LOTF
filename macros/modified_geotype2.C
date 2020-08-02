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
  PndGeoHandling geoH(1572086365,"evtcomplete_par.root");
  geoH.FillSensorMap();

  //FairBaseParSet *parbase = f.Get("FairBaseParSet");

  PndGeoSttPar *parameters =  (PndGeoSttPar*) f.Get("PndGeoSttPar");
  int fGeoType = parameters->GetGeometryType();
   cout << fGeoType << endl;
// 

  PndSttMapCreator *mapper = new PndSttMapCreator(parameters);
 parameters->SetGeometryType(1);
cout << "HERE"<<endl;
  TClonesArray *fTubeArray = mapper->FillTubeArray();
PndSttGeometryMap *fMap; //!

  map<int, vector< vector<int> > > fStrawIndex;
  map<int, TArrayI> fStrawNeighbors;
  vector<int> fSectorOfStraw;
  vector<int> fLayerOfStraw;
  vector<bool> fAxialStraw;
  vector<double> fSectorStart;
  vector<double> fSectorEnd;
  bool fStrawMapInitialized;

cout << "Does it work ?"<<endl;
fSectorStart.push_back(1.57);
fSectorEnd.push_back(2.62);
fSectorStart.push_back(2.62);
fSectorEnd.push_back(-2.62+2*TMath::Pi());
fSectorStart.push_back(-2.62+2*TMath::Pi());
fSectorEnd.push_back(-1.57+2*TMath::Pi());
fSectorStart.push_back(-1.57+2*TMath::Pi());
fSectorEnd.push_back(-0.52+2*TMath::Pi());
fSectorStart.push_back(-0.52+2*TMath::Pi());
fSectorEnd.push_back(0.52);
fSectorStart.push_back(0.52);
fSectorEnd.push_back(1.57);
 vector<int> currentRow;
  int sector = 0;
  int lastsector = 0;
  int row = 0;
int fVerbose = 1;
  fSectorOfStraw.push_back(-1);
  fLayerOfStraw.push_back(-1);
  fAxialStraw.push_back(false);
PndSttTube* tube = NULL;
 cout << "Generating straw map for " << fTubeArray->GetEntriesFast() << " straws." << endl;
  for (int i = 1; i < fTubeArray->GetEntriesFast(); i++) {
   cout << "Finding tube." << endl;
   tube = (PndSttTube*) fTubeArray->At(i);
     cout << "Tube address: " << tube << endl;
    bool isaxial = (tube->GetWireDirection().Theta() < 0.001);

      cout << "Axial Straw: " << isaxial << endl;
      cout << "Wire Direction: ";
      tube->GetWireDirection().Print();
    
    double phi = tube->GetPosition().Phi();
    cout << "Phi: " << phi << endl;
    if (phi < 0) phi += 2*TMath::Pi();
     cout << "Checking sector." << endl;
    while ( !( (phi > fSectorStart.at(sector)) && (phi < fSectorEnd.at(sector)) )
	    && !( (sector == 4) && ( (phi > fSectorStart.at(4)) || (phi < fSectorEnd.at(4)) ) ) ) {
      sector++;
      sector %= 6;
      }
    if (sector != lastsector) {
      fStrawIndex[lastsector].push_back(currentRow);
      currentRow.clear();
      if(fVerbose > 0) cout << "SECTOR COMPLETE: Row " << row << " added to sector " << lastsector << endl;
    }
    if (sector < lastsector) row++;
    lastsector = sector;
    currentRow.push_back(i);
    fSectorOfStraw.push_back(sector);
    fLayerOfStraw.push_back(row);
    fAxialStraw.push_back(isaxial);
    if(fVerbose > 0) cout << "Straw " << i << " added to " << sector << ", " << row << endl;

    tube->SetLayerID(row);
    tube->SetSectorID(sector); 
    

  
  }


// TCanvas *c = new TCanvas("c", "", 400, 0, 900, 900);
// TH2F *h2 = new TH2F("h2", "", 100, -42, 42, 100, -42, 42);
//  h2->SetStats(kFALSE);
//  h2->Draw();
  
  TArc *arc = NULL;
  TMarker *mrk = NULL;
  //  PndSttTube *tube = NULL;
PndSttTube *sttTube = nullptr;
int nodeID;
int nodeType;
int SectorLimit = 0;
 int Sector = 0;
   int Layer = 0;
  bool LayerLimit = 0 ;
  int itube = 0;

for(itube = 1; itube < fTubeArray->GetEntriesFast(); itube++) {
  tube = (PndSttTube*) fTubeArray->At(itube);
  //sttTube = PndSttTubeMap::Instance()->GetTube(iTube);

  /* if(tube->IsParallel()) {
    arc = new TArc(tube->GetPosition().X(), tube->GetPosition().Y(), 0.5);
    arc->Draw("SAME");
  }
  else {
    mrk = new TMarker(tube->GetPosition().X(), tube->GetPosition().Y(), 6);
    mrk->Draw("SAME");
  }
  */
  Sector     = tube->GetSectorID();
  Layer      = tube->GetLayerID();
  LayerLimit = tube->IsLayerLimit();
  SectorLimit = tube->IsSectorLimit();
  cout << tube->IsLayerLimit() << "\t";
 }
/* while(itube != -1) {
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
      cout << Sector << endl;
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
  }*/
}
