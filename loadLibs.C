/* ********************************
 * Author: M. Babai               *
 * Version:                       *
 * Licence:                       *
 * ********************************/

Bool_t isLibrary(const char* libName)
{
  if (TString(gSystem->DynamicPathName(libName, kTRUE)) != TString(""))
    return kTRUE;
  else  
    return kFALSE;
}

//rootlogon()
void loadLibs()
{
  // gROOT->LoadMacro("$VMCWORKDIR/gconfig/basiclibs.C");
  // basiclibs();
  gSystem->AddIncludePath("-I$VMCWORKDIR/fairtools/");
  
  gSystem->AddIncludePath("-I$VMCWORKDIR/base/");
  gSystem->AddIncludePath("-I$VMCWORKDIR/base/sim/");
  gSystem->AddIncludePath("-I$VMCWORKDIR/base/steer/");
  gSystem->AddIncludePath("-I$VMCWORKDIR/base/event/");
  
  gSystem->AddIncludePath("-I$VMCWORKDIR/dbase/");
  gSystem->AddIncludePath("-I$VMCWORKDIR/dbase/dbInterface/");
  gSystem->AddIncludePath("-I$VMCWORKDIR/dbase/dbInput/");
  gSystem->AddIncludePath("-I$VMCWORKDIR/dbase/dbUtils/");
  gSystem->AddIncludePath("-I$VMCWORKDIR/dbase/dbValidation/");
  
  gSystem->AddIncludePath("-I$VMCWORKDIR/emc/EmcData/");
  gSystem->AddIncludePath("-I$VMCWORKDIR/emc/EmcCorr/");
  gSystem->AddIncludePath("-I$VMCWORKDIR/emc/EmcDigi/");
  gSystem->AddIncludePath("-I$VMCWORKDIR/emc/EmcMC/");
  gSystem->AddIncludePath("-I$VMCWORKDIR/emc/EmcReco/");
  gSystem->AddIncludePath("-I$VMCWORKDIR/emc/EmcTools/");
  
  gSystem->AddIncludePath("-I$VMCWORKDIR/pndbase/PndStdUtils/");
  gSystem->AddIncludePath("-I$VMCWORKDIR/pnddata/");
  gSystem->AddIncludePath("-I$VMCWORKDIR/pnddata/TrackData/");
  gSystem->AddIncludePath("-I$VMCWORKDIR/pnddata/SdsData/");
  
  gSystem->AddIncludePath("-I$VMCWORKDIR/trackbase/");
  gSystem->AddIncludePath("-I$VMCWORKDIR/tracking/");
  
  gSystem->AddIncludePath("-I$VMCWORKDIR/parbase/");
  gSystem->AddIncludePath("-I$VMCWORKDIR/stt/");
  gSystem->AddIncludePath("-I$VMCWORKDIR/pnddata/SttData/");
  gSystem->AddIncludePath("-I$VMCWORKDIR/stt/sttreco");
  
  // Load Panda libraries
  if(isLibrary("libFairTools"))gSystem->Load("libFairTools");
  if(isLibrary("libFairDB")) gSystem->Load("libFairDB");
  if(isLibrary("libGeoBase"))gSystem->Load("libGeoBase");
  if(isLibrary("libParBase"))gSystem->Load("libParBase");
  if(isLibrary("libBase"))   gSystem->Load("libBase");
  if(isLibrary("libMCStack")) gSystem->Load("libMCStack");
  
  if(isLibrary("libDpmEvtGen"))gSystem->Load("libDpmEvtGen");
  if(isLibrary("libpythia8"))gSystem->Load("libpythia8");
  if(isLibrary("libFlukaResults"))gSystem->Load("libFlukaResults");
  if(isLibrary("libPhotos")){// these three depend on each other
    gSystem->Load("libPhotos");
    if(isLibrary("libEvtGen")){
      gSystem->Load("libEvtGen");
      if(isLibrary("libEvtGenDirect"))gSystem->Load("libEvtGenDirect");
    }
  }
  if(isLibrary("libPndBase"))gSystem->Load("libPndBase");
  if(isLibrary("libGlobalTasks"))gSystem->Load("libGlobalTasks");
  if(isLibrary("libTrkBase"))gSystem->Load("libTrkBase");
  if(isLibrary("libgeneralTools"))gSystem->Load("libgeneralTools");
  if(isLibrary("libPndData"))gSystem->Load("libPndData");
  if(isLibrary("libbuffers"))gSystem->Load("libbuffers");
  if(isLibrary("libField"))gSystem->Load("libField");
  if(isLibrary("libPassive"))gSystem->Load("libPassive");
  if(isLibrary("libGen"))gSystem->Load("libGen");
  if(isLibrary("libPGen"))gSystem->Load("libPGen");
  if(isLibrary("libEmc"))gSystem->Load("libEmc"); 
  if(isLibrary("libgenfit"))gSystem->Load("libgenfit");
  if(isLibrary("libtrackrep"))gSystem->Load("libtrackrep");
  if(isLibrary("libgenfitAdapters"))gSystem->Load("libgenfitAdapters");
  if(isLibrary("libriemann"))gSystem->Load("libriemann");
  if(isLibrary("libTpcBase"))gSystem->Load("libTpcBase"); 
  if(isLibrary("libTpc"))gSystem->Load("libTpc"); 
  if(isLibrary("libTpcReco"))gSystem->Load("libTpcReco");
  if(isLibrary("libStt"))gSystem->Load("libStt");
  if(isLibrary("libSttReco"))gSystem->Load("libSttReco");
  if(isLibrary("libSds"))gSystem->Load("libSds");
  if(isLibrary("libSdsReco"))gSystem->Load("libSdsReco");
  if(isLibrary("libMvd"))gSystem->Load("libMvd");
  if(isLibrary("libMvdReco"))gSystem->Load("libMvdReco");
  if(isLibrary("libMvdTrk"))gSystem->Load("libMvdTrk");
  if(isLibrary("libSttMvdTracking"))gSystem->Load("libSttMvdTracking");
  if(isLibrary("libGem"))gSystem->Load("libGem");
  if(isLibrary("libTof"))gSystem->Load("libTof");
  if(isLibrary("libDrcProp"))gSystem->Load("libDrcProp");
  if(isLibrary("libDrc"))gSystem->Load("libDrc");
  if(isLibrary("libMdt"))gSystem->Load("libMdt");
  if(isLibrary("libDch"))gSystem->Load("libDch");
  if(isLibrary("libLheTrack"))gSystem->Load("libLheTrack");
  if(isLibrary("libGeane"))gSystem->Load("libGeane");
  if(isLibrary("libRpc"))gSystem->Load("libRpc");
  if(isLibrary("libLumi"))gSystem->Load("libLumi");
  if(isLibrary("libRho"))gSystem->Load("libRho");
  if(isLibrary("libTMVA"))gSystem->Load("libTMVA.so");
  if(isLibrary("libAnalysisTools"))gSystem->Load("libAnalysisTools"); 
  if(isLibrary("libPid"))gSystem->Load("libPid");
  if(isLibrary("librecotasks"))gSystem->Load("librecotasks");
  if(isLibrary("libRecoHits"))gSystem->Load("libRecoHits");
  if(isLibrary("libRecoTasks"))gSystem->Load("libRecoTasks");
  if(isLibrary("libEnDrc"))gSystem->Load("libEnDrc");
  if(isLibrary("libDsk"))gSystem->Load("libDsk");
  if(isLibrary("libGlobal"))gSystem->Load("libGlobal");
  if(isLibrary("libMCMatch"))gSystem->Load("libMCMatch");
  if(isLibrary("libMva"))gSystem->Load("libMva");

  //gSystem->Load("libMemStat");
}


//////////////////// EXTRA DOCUMENTATION ////////////////////////
// .L MyScript.C+g
// .L MyScript.C+O
// gSystem->SetAclicMode(TSystem::kDebug);
// gSystem->SetAclicMode(TSystem::kOpt);
// gSystem->AddIncludePath(" -I$HOME/mypackage/include ")
// You can also overwrite the existing include path:
// gSystem->SetIncludePath(" -I$HOME/mypackage/include ")
// To add library that should be used during linking of the shared library use something like:
// gSystem->AddLinkedLibs("-L/my/path -lanylib");

//For example, in one script you can use ACLiC to compile and load another script.
//gROOT->ProcessLine(".L MyScript.C+")
//gROOT->ProcessLine(".L MyScript.C++")
// OR we can use
//gROOT->LoadMacro("pathopen.cpp+");
//////////////////////////////////////////////////////
