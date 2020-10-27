CXX=g++
OBJDIR=build
DATADIR=data
BINDIR=bin
SRCDIR=src
INCDIR=include

ROOTFLAGS=$(shell ${ROOTSYS}/bin/root-config --cflags)
ROOTLIBS=$(shell ${ROOTSYS}/bin/root-config --libs)
CXXFLAGS=-std=c++11 -Wall -O3 -I./$(INCDIR)

.PHONY: all
all: $(BINDIR)/tcs-gen $(BINDIR)/params-inspect $(BINDIR)/root-to-dat $(BINDIR)/brem-compare $(BINDIR)/integrate $(BINDIR)/GenOptions.dat $(BINDIR)/CFFs_DD_Feb2012.dat

.PHONY: clean
clean:
	rm -rf $(OBJDIR)/*
	rm -rf $(BINDIR)/*

$(OBJDIR):
	mkdir -p $(OBJDIR)
$(BINDIR):
	mkdir -p $(BINDIR)

$(BINDIR)/tcs-gen: $(SRCDIR)/TCSGen.cc $(OBJDIR)/TTCSKine.o $(OBJDIR)/KinFunctions.o $(OBJDIR)/CrsFunctions.o $(OBJDIR)/TTCSCrs.o $(OBJDIR)/GPDs.o | $(BINDIR)
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(ROOTLIBS) -o $(BINDIR)/tcs-gen $(SRCDIR)/TCSGen.cc $(OBJDIR)/TTCSKine.o $(OBJDIR)/KinFunctions.o $(OBJDIR)/CrsFunctions.o $(OBJDIR)/TTCSCrs.o $(OBJDIR)/GPDs.o

$(BINDIR)/params-inspect: $(SRCDIR)/ParamsInspect.cc | $(BINDIR)
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(ROOTLIBS) -o $(BINDIR)/params-inspect $(SRCDIR)/ParamsInspect.cc

$(BINDIR)/root-to-dat: $(SRCDIR)/RootToDat.cc | $(BINDIR)
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(ROOTLIBS) -o $(BINDIR)/root-to-dat $(SRCDIR)/RootToDat.cc

$(BINDIR)/brem-compare: $(SRCDIR)/BremCompare.cc $(OBJDIR)/KinFunctions.o | $(BINDIR)
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(ROOTLIBS) -o $(BINDIR)/brem-compare $(SRCDIR)/BremCompare.cc $(OBJDIR)/KinFunctions.o

$(BINDIR)/integrate: $(SRCDIR)/Integrate.cc | $(BINDIR)
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(ROOTLIBS) -o $(BINDIR)/integrate $(SRCDIR)/Integrate.cc

$(BINDIR)/CFFs_DD_Feb2012.dat: $(DATADIR)/CFFs_DD_Feb2012.dat | $(BINDIR)
	cp $(DATADIR)/CFFs_DD_Feb2012.dat $(BINDIR)/CFFs_DD_Feb2012.dat

$(BINDIR)/GenOptions.dat: $(DATADIR)/GenOptionsStandard.dat | $(BINDIR)
	cp $(DATADIR)/GenOptionsStandard.dat $(BINDIR)/GenOptions.dat

$(OBJDIR)/TTCSKine.o: $(SRCDIR)/TTCSKine.cc $(INCDIR)/TTCSKine.h | $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c -o $(OBJDIR)/TTCSKine.o $(SRCDIR)/TTCSKine.cc

$(OBJDIR)/GPDs.o: $(SRCDIR)/GPDs.cc $(INCDIR)/GPDs.h | $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c -o $(OBJDIR)/GPDs.o $(SRCDIR)/GPDs.cc

$(OBJDIR)/TTCSCrs.o: $(SRCDIR)/TTCSCrs.cc $(INCDIR)/TTCSCrs.h $(OBJDIR)/GPDs.o | $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c -o $(OBJDIR)/TTCSCrs.o $(SRCDIR)/TTCSCrs.cc $(OBJDIR)/GPDs.o

$(OBJDIR)/KinFunctions.o: $(SRCDIR)/KinFunctions.cc $(INCDIR)/KinFunctions.h | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c -o $(OBJDIR)/KinFunctions.o $(SRCDIR)/KinFunctions.cc

$(OBJDIR)/CrsFunctions.o: $(SRCDIR)/CrsFunctions.cc $(INCDIR)/CrsFunctions.h | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c -o $(OBJDIR)/CrsFunctions.o $(SRCDIR)/CrsFunctions.cc

