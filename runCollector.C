/* ********************************
 * Author: M. Babai               *
 * Version:                       *
 * Licence:                       *
 * ********************************/
#include "loadLibs.C"

void runCollector( std::string const &OutFileName = "LOTFAnalysis.root",		   
		   int firstEvt =0, //9, //0,190
		   int lastEvt =10
		   )//50->61, 4->5
{
  // Load basic libs and headers.
  // gROOT->LoadMacro("loadLibs.C");
  loadLibs();

  // Compile and load nedded modules.
  gSystem->CompileMacro("auxfunctions.cpp","kO");
  gSystem->CompileMacro("gridNode.cpp","kO");
  gSystem->CompileMacro("hitcoordinate.cpp","kO");
  gSystem->CompileMacro("CoordGrid.cpp","kO");
  gSystem->CompileMacro("pathCandidate.cpp","kO");
  gSystem->CompileMacro("trackObject.cpp","kO");
  gSystem->CompileMacro("SttMVDEventDataReader.cpp","kO");
  gSystem->CompileMacro("logc.cpp","kO");
  gSystem->CompileMacro("CollectSttMvdPoints.cpp","kO");
  gSystem->CompileMacro("phconnect.cpp","kO");
  gSystem->CompileMacro("phfitting.cpp","kO");
  gSystem->CompileMacro("phmerging.cpp","kO");
  gSystem->CompileMacro("phzinterp.cpp","kO");
  gSystem->CompileMacro("floodingFilter.cpp","kO");
  gSystem->CompileMacro("error.cpp","kO");
  gSystem->CompileMacro("circle.cpp","kO");
  gSystem->CompileMacro("CreatePlotAllEventComponents.C","kO");

  gROOT->ProcessLine(Form("floodingFilter(\"%s\", %d, %d);", OutFileName.c_str(), firstEvt, lastEvt));


  // OLD CODE
    //  gSystem->CompileMacro("mvdMapCreator.cpp","kO");
  // gSystem->CompileMacro("auxiliaryfunctions.cpp","kO");
  //gSystem->CompileMacro("utilfunctions.cpp","kO");
  // gSystem->CompileMacro("path_queue.cpp","kO");
  //  gSystem->CompileMacro("pathopen.cpp","kO");
  //  gSystem->CompileMacro("queue.cpp","kO");
  // gSystem->CompileMacro("reconstruction.cpp","kO");

  //======= Load the Plot and helper macro's ====
  //   gSystem->CompileMacro("PlotTracks.C","k0");
   //  gSystem->CompileMacro("PlotGrid.C","k0");
  //  gSystem->CompileMacro("PlotOrientations.C","kO");
  //  gSystem->CompileMacro("PlotConnectedComponents.C","kO");
  //  gSystem->CompileMacro("PlotMergedComponents.C","kO");
  //gSystem->CompileMacro("SttHitCheck.C","k0");
  // gSystem->CompileMacro("DrawCorrected.C","k0");
  // gSystem->CompileMacro("ComputePerTrackError.C","k0");
  //==================
  // SttHitCheck("sttHitControleHist.root","DummyOut.root");
  /* Collect event data and geometry parameters from simulated data*/
  //  gSystem->CompileMacro("SttMVDPointsToTreeFileWriter.cpp","kO");
  //  SttMVDPointsToTreeFileWriter();
  //gROOT->ProcessLine(Form("TreeToOut(\"%s\", \"%s\",\"%s\",\"%s\" ,\"%s\", \"%s\");",SimInFile.c_str(),digiInFile.c_str(),SimParamInfile.c_str(),OutFileName.c_str(),TreeName.c_str(),nodetree.c_str()));
  /* Perform the actual filtering (run main) */
  // gROOT->ProcessLine(Form("performFilter(%ld, %ld, %ld, %f, %f, %ld, \"%s\", %d, %d);", plength, area, MinResponce, lambda, tol, gapSize,OutFileName.c_str(), firstEvt, lastEvt));

  // gROOT->ProcessLine(Form("CreatePlotAllEventComponents(\"%s\");", OutFileName.c_str()));

  //_________ Make plots ________
  
  /* PlotGrid( Infile, ww, hh, dim )*/
  //PlotGrid(OutFileName, WindowSize, WindowSize, 2);

  /*PlotTracks(Infile, ww, hh,nEv, dim)*/
  //_____ PlotTracks(OutFileName, WindowSize, WindowSize, numEvt, dim);

  /*PlotOrientations(Infile, ww, hh, nEv, dim, nOrient);*/
  //PlotOrientations(OutFileName, WindowSize, WindowSize, numEvt, dim, Orientations);
  
  /* PlotConnectedComponents(Infile, ww, hh, nEv, nComp, dim, showGrid);*/
  // PlotConnectedComponents(OutFileName, WindowSize, WindowSize, numEvt, numConCom, dim, true);
  // PlotConnectedComponents(OutFileName, WindowSize, WindowSize, numEvt, numConCom, 2, false);
  
  // MVD Merged components
  // PlotMergedComponents(OutFileName, WindowSize, WindowSize, numEvt, numConCom, dim, true);

  // PlotMergedComponents(OutFileName, WindowSize, WindowSize, numEvt, numConCom, dim, false);

  // Write event Plots
  
  /*  TFile outFile(OutFileName.c_str(), "UPDATE");
  outFile.mkdir("CanvasPlots");
  outFile.Close();
  //____ Make plots for all events and connected components
  
  if( (firstEvt < 0) && (lastEvt < 0) ) {
    firstEvt = 0;
    lastEvt  = numEvt;
  }
  for(size_t e = 0; e < (lastEvt - firstEvt); ++e) {
    gROOT->ProcessLine(Form("CreatePlotAllEventComponents(\" %s\", \"CanvasPlots\", %ld, %ld, %ld, %ld", OutFileName.c_str(), WindowSize, WindowSize, dim, e));
    }*/
  // Plot Corrected coordinates for skewed detectors.
  //DrawCorrected();
}
//________________________________________________________
//gROOT->LoadMacro("TestSttCoord.C+");
// TestSttCoord();
//gROOT->ProcessLine(".L gridNode.cpp+O");
// gROOT->ProcessLine(".L hitcoordinate.cpp+O");
// gROOT->ProcessLine(".L CoordGrid.cpp+O");
// gROOT->ProcessLine(".L utilfunctions.cpp+O");
// gROOT->ProcessLine(".L CollectSttMvdPoints.cpp+O");
// gROOT->ProcessLine(".L path_queue.cpp+O");
// gROOT->ProcessLine(".L trackObject.cpp+O");
// gROOT->ProcessLine(".L pathopen.cpp+O");
// gROOT->ProcessLine(".L performFilter.cpp+O");
  
//======= Load the Plot macro's ====
// gROOT->LoadMacro("PlotTracks.C");
// gROOT->LoadMacro("PlotGrid.C");
// gROOT->ProcessLine(".L PlotOrientations.C+O");
// gROOT->ProcessLine(".L PlotConnectedComponents.C+O");
// gROOT->ProcessLine(".L PlotMergedComponents.C+O");
/*
gSystem->CompileMacro(fileName,options);
The possible options are:
k : keep the shared library after the session end.
f : force recompilation.
g : compile with debug symbol
O : optimized the code (ignore if 'g' is specified)
c : compile only, do not attempt to load the library.
- : if buildir is set, use a flat structure (see buildir below)
*/
