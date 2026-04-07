# Simple build/run helper for Tutorial_AI
#
# Notes:
# - AIZ.C is a ROOT macro (no main()), so we do not build it as a standalone executable.
# - We build TLVUtils as an object and shared library, then run AIZ via ROOT.

CXX        ?= g++
ROOTCONFIG ?= root-config

ROOTFLAGS := $(shell $(ROOTCONFIG) --cflags)
ROOTLIBS  := $(shell $(ROOTCONFIG) --libs)
CXXFLAGS  ?= -O2 -g -fPIC

TARGET_LIB = libTLVUtils.so
OBJECTS    = TLVUtils.o

.PHONY: all libs run run-y run-test clean

all: libs

libs: $(TARGET_LIB)

TLVUtils.o: TLVUtils.cxx TLVUtils.h
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

$(TARGET_LIB): TLVUtils.o
	$(CXX) -shared -o $@ $^ $(ROOTLIBS)

# Run AIZ in pT mode (AIZ(false))
run: libs
	root -l -b -q -e '.L TLVUtils.cxx' -e '.L AIZ.C' -e 'AIZ(false);' -e '.q'

# Run AIZ in |y| mode (AIZ(true))
run-y: libs
	root -l -b -q -e '.L TLVUtils.cxx' -e '.L AIZ.C' -e 'AIZ(true);' -e '.q'

# Quick test mode
run-test: libs
	root -l -b -q -e '.L TLVUtils.cxx' -e '.L AIZ.C' -e 'test=true;' -e 'AIZ(false);' -e '.q'

clean:
	rm -f $(OBJECTS) $(TARGET_LIB)
