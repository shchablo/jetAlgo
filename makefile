in=main

out=$(in)
bin_name=bin/$(out)
source_name=$(in).C

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --glibs)

FASTJETFLAGS := $(shell fastjet-install/bin/fastjet-config --cxxflags)
FASTJETLIBS  := $(shell fastjet-install/bin/fastjet-config --libs --plugins)

all: 

	g++ $(source_name) source/*.cc -I include/ $(ROOTFLAGS) $(ROOTLIBS) $(FASTJETFLAGS) $(FASTJETLIBS) -o jetAlgo
