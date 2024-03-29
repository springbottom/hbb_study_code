# --------------------------------------------- #
# Makefile for myexample code                   #
# Pascal Nef, updated by Ben Nachman (2017)     #
#                                               #
# Note: source setup.sh before make             #
# --------------------------------------------- #

CXXFLAGS = -g -std=c++0x  -O2 -Wall 
LDFLAGS = -g

.PHONY: clean debug all

all: myexample

myexample:  myexample.so  myexampleAnalysis.so 
	$(CXX) myexample.so myexampleAnalysis.so -o $@.exe \
	$(CXXFLAGS) -Wno-shadow  \
	`root-config --glibs`  \
	-L$(FASTJETLOCATION)/lib `$(FASTJETLOCATION)/bin/fastjet-config --libs --plugins ` -lRecursiveTools \
	-L$(PYTHIA8LOCATION)/lib -lpythia8 -ldl \
	-L$(BOOSTLIBLOCATION) -lboost_program_options 

myexample.so: myexample.C  myexampleAnalysis.so   
	$(CXX) -o $@ -c $<  \
	$(CXXFLAGS) -Wno-shadow -fPIC -shared \
	`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins` \
	-I$(PYTHIA8LOCATION)/include \
	-I $(BOOSTINCDIR) \
	`root-config --cflags` 

myexampleAnalysis.so : myexampleAnalysis.cc myexampleAnalysis.h
	$(CXX) -o $@ -c $<  \
	$(CXXFLAGS) -Wno-shadow -fPIC -shared \
	`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins`  \
	-I$(PYTHIA8LOCATION)/include \
	`root-config --cflags`   

clean:
	rm -rf *.exe
	rm -rf *.o
	rm -rf *.so
	rm -f *~
