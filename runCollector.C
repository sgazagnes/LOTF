/* ********************************
 * Author: M. Babai               *
 * Version:                       *
 * Licence:                       *
 * ********************************/
#include "loadLibs.C"

void runCollector( std::string const &OutFileName = "LOTFAnalysis.root",		   
		   int firstEvt =0, //9, //0,190
		   int lastEvt =20000
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
  gSystem->CompileMacro("error.cpp","kO");
  gSystem->CompileMacro("circle.cpp","kO");
  gSystem->CompileMacro("CreatePlotAllEventComponents.C","kO");
  
  gSystem->CompileMacro("floodingFilter.cpp","kO");

  gROOT->ProcessLine(Form("floodingFilter(\"%s\", %d, %d);", OutFileName.c_str(), firstEvt, lastEvt));
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
