#/*************************************
# * Author: M. Babai (M.Babai@rug.nl) *
# * Version:                          *
# * License:                          *
# HOW TO RUN                          *
# OMP_NUM_THREADS=4                   *
#                                     *
#/*************************************
# Setting the linker and compiler flags and options
#####################################################
FLAGSMACHINE = -march=native -mtune=native

CCFLAGS  = -W -Wall -Wextra -Wredundant-decls -Wshadow #-Werror
CCFLAGS += -Woverlength-strings  -Wfloat-equal #-Weffc++
CCFLAGS += -O2 -fpic -ansi -pedantic -std=c++11
CCFLAGS += -fPIE #-fstack-check -funroll-loops
#CCFLAGS += -fopenmp -malign-double
CCFLAGS += ${FLAGSMACHINE}
CXX = g++
INCLUDE = $(shell root-config --cflags) -I./

LIBS  = $(shell root-config --libs) -lTreePlayer -lXMLIO -fopenmp
# LIBS += -lTMVA -lMinuit -lRooFitCore -lFoam
# LIBS += -lMLP

TRGTS  =  gridNode.o pathCandidate.o hitcoordinate.o CoordGrid.o utilfunctions.o
TRGTS  += path_queue.o trackObject.o pathopen.o SttMVDEventDataReader.o
TRGTS  += performFilter.o auxiliaryfunctions.o runTracking.o executable

all: $(TRGTS)

auxiliaryfunctions.o: auxiliaryfunctions.cpp
	$(CXX) $(CCFLAGS) $(INCLUDE) -c $< -o $@

gridNode.o: gridNode.cpp
	$(CXX) $(CCFLAGS) $(INCLUDE) -c $< -o $@

pathCandidate.o: pathCandidate.cpp
	$(CXX) $(CCFLAGS) $(INCLUDE) -c $< -o $@

hitcoordinate.o: hitcoordinate.cpp
	$(CXX) $(CCFLAGS) $(INCLUDE) -c $< -o $@

CoordGrid.o: CoordGrid.cpp
	$(CXX) $(CCFLAGS) $(INCLUDE) -c $< -o $@

utilfunctions.o: utilfunctions.cpp
	$(CXX) $(CCFLAGS) $(INCLUDE) -c $< -o $@

path_queue.o: path_queue.cpp
	$(CXX) $(CCFLAGS) $(INCLUDE) -c $< -o $@

trackObject.o: trackObject.cpp
	$(CXX) $(CCFLAGS) $(INCLUDE) -c $< -o $@

pathopen.o: pathopen.cpp
	$(CXX) $(CCFLAGS) $(INCLUDE) -c $< -o $@

SttMVDEventDataReader.o: SttMVDEventDataReader.cpp
	$(CXX) $(CCFLAGS) $(INCLUDE) -c $< -o $@

performFilter.o: performFilter.cpp
	$(CXX) $(CCFLAGS) $(INCLUDE) -c $< -o $@

runTracking.o: runTracking.cpp
	$(CXX) $(CCFLAGS) $(INCLUDE) -c $< -o $@

executable: gridNode.o pathCandidate.o hitcoordinate.o CoordGrid.o utilfunctions.o \
	path_queue.o trackObject.o pathopen.o SttMVDEventDataReader.o \
	performFilter.o auxiliaryfunctions.o runTracking.o
	$(CXX) $(CCFLAGS) $(LIBS) trackObject.o path_queue.o hitcoordinate.o \
	gridNode.o utilfunctions.o pathCandidate.o CoordGrid.o pathopen.o \
	SttMVDEventDataReader.o auxiliaryfunctions.o \
	runTracking.o -o runTracking

## Cleaning and remove targets
clean:
	rm -rf *~ *.d *.o
	$(MAKE) -C macros clean;

clean_objects: clean
	rm -rf *.so *.o runTracking

clean_root:
	rm -rf *.root

clean_all: clean clean_root
	rm -rf *.so *.log *.o runTracking

distclean: clean_all
	rm -rf *.root *.pdf *.eps exampleCoordinates_CSV.txt \
	GridDumpFile.txt runTracking
	$(MAKE) -C macros distclean;
