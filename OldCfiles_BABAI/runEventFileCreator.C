/* ********************************
 * Author: M. Babai               *
 * Version:                       *
 * Licence:                       *
 * ********************************/
#include "loadLibs.C"

void runEventFileCreator(
			 std::string const &SimInFile= "rootfiles/evtcomplete_sim.root",// Sim input file
			 std::string const &digiInFile= "rootfiles/evtcomplete_digi.root",// Digi input file
			 std::string const &SimParamInfile= "rootfiles/evtcomplete_par.root",// Simparam inputfile
			 std::string const &OutFileName= "rootfiles/CompactEventData.root",// Outputfile name
			 std::string const &TreeName= "EventData", // Output tree name
			 std::string const &nodetree= "Geometry"
			  )
{
  // Load basic libs and headers.
  //gROOT->LoadMacro("loadLibs.C");
  loadLibs();
  
  // define the Compiler flags
  std::string CompileOptions = "-O3 -march=native -mtune=native -fPIC";
  // CompileOptions += " -funroll-loops -fopenmp -malign-double";
  // CompileOptions += " -Wshadow -Weffc++ -W -Wall -Wextra -Wfloat-equal";
  // CompileOptions += " -Wredundant-decls";
  
  gSystem->SetAclicMode(TSystem::kOpt);
  gSystem->SetFlagsOpt(CompileOptions.c_str());

  // Compile and load nedded modules.
  gSystem->CompileMacro("SimToCompactTree.cpp","kO");
  gROOT->ProcessLine(Form("TreeToOut(\"%s\", \"%s\",\"%s\",\"%s\" ,\"%s\", \"%s\");",SimInFile.c_str(),digiInFile.c_str(),SimParamInfile.c_str(),OutFileName.c_str(),TreeName.c_str(),nodetree.c_str()));

  exit(0);
}
