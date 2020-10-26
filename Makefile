CXX=g++
OBJDIR=build
BINDIR=bin
SRCDIR=src
INCDIR=include

ROOTFLAGS=$(shell ${ROOTSYS}/bin/root-config --cflags)
ROOTLIBS=$(shell ${ROOTSYS}/bin/root-config --libs)
CXXFLAGS=-std=c++11 -Wall -I./$(INCDIR)

.PHONY: clean
clean:
	rm -rf $(OBJDIR)/*
	rm -rf $(BINDIR)/*

.PHONY: all
all: $(BINDIR)/tcs-gen $(BINDIR)/root-to-dat

$(BINDIR)/tcs-gen: $(SRCDIR)/TCSGen.cc $(OBJDIR)/TTCSKine.o $(OBJDIR)/KinFunctions.o $(OBJDIR)/CrsFunctions.o $(OBJDIR)/TTCSCrs.o $(OBJDIR)/GPDs.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(ROOTLIBS) -o $(BINDIR)/tcs-gen $(SRCDIR)/TCSGen.cc $(OBJDIR)/TTCSKine.o $(OBJDIR)/KinFunctions.o $(OBJDIR)/CrsFunctions.o $(OBJDIR)/TTCSCrs.o $(OBJDIR)/GPDs.o

$(BINDIR)/root-to-dat: $(SRCDIR)/RootToDat.cc
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(ROOTLIBS) -o $(BINDIR)/root-to-dat $(SRCDIR)/RootToDat.cc

$(OBJDIR)/TTCSKine.o: $(SRCDIR)/TTCSKine.cc $(INCDIR)/TTCSKine.h
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c -o $(OBJDIR)/TTCSKine.o $(SRCDIR)/TTCSKine.cc

$(OBJDIR)/GPDs.o: $(SRCDIR)/GPDs.cc $(INCDIR)/GPDs.h
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c -o $(OBJDIR)/GPDs.o $(SRCDIR)/GPDs.cc

$(OBJDIR)/TTCSCrs.o: $(SRCDIR)/TTCSCrs.cc $(INCDIR)/TTCSCrs.h $(OBJDIR)/GPDs.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c -o $(OBJDIR)/TTCSCrs.o $(SRCDIR)/TTCSCrs.cc $(OBJDIR)/GPDs.o

$(OBJDIR)/KinFunctions.o: $(SRCDIR)/KinFunctions.cc $(INCDIR)/KinFunctions.h
	$(CXX) $(CXXFLAGS) -c -o $(OBJDIR)/KinFunctions.o $(SRCDIR)/KinFunctions.cc

$(OBJDIR)/CrsFunctions.o: $(SRCDIR)/CrsFunctions.cc $(INCDIR)/CrsFunctions.h
	$(CXX) $(CXXFLAGS) -c -o $(OBJDIR)/CrsFunctions.o $(SRCDIR)/CrsFunctions.cc

