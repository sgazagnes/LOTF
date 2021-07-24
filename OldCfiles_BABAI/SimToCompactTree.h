/*************************************
 * Author: M. Babai (M.Babai@rug.nl) *
 * Version:                          *
 * License:                          *
 *************************************/
#pragma once
#ifndef SIMTOCOMPACTTREE_H
#define SIMTOCOMPACTTREE_H

// C and C++ standard headers
#include <string>

// Root headers
double Cartesian_To_Polar(float const x, float const y,
                          std::pair<float,float>& polarOut,
                          bool useSign=true);

void TreeToOut(std::string const &SimInFile,// Sim input file
				     std::string const &digiInFile,// Digi input file
				     std::string const &SimParamInfile,// Simparam inputfile
				     std::string const &OutFileName,// Outputfile name
				     std::string const &TreeNameEvt, // Output tree name
				     std::string const &NodeTreeName // Tree name to hold node name
				     );
#endif
