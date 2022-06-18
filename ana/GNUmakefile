OBJS = meTree.o analysis.o
DEPS = meTree.h
EXE = analysis

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)

#INCFLAGS = -I$(ROOTSYS)/include
#LDFLAGS = -L$(ROOTSYS)/lib
ROOTINC = $(shell root-config --incdir)

CXX = g++
CFLAGS = -m32 -Wall -g $(ROOTCFLAGS)

$(EXE): $(OBJS)
	$(CXX) -o $(EXE) $(OBJS) $(ROOTCFLAGS) $(ROOTLIBS)

%.o: %.cxx
	$(CXX) $(CFLAGS) -c -o $@  $<

.PHONY: clean
clean:
	rm -f ${EXE} ${OBJS}
