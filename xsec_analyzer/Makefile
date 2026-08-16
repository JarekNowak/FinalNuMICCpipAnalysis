# Verify that ROOT is set up before attempting a build. Don't bother with the
# check if we just want to "make clean"
ifneq ($(MAKECMDGOALS),clean)
  ROOTCONFIG := $(shell command -v root-config 2> /dev/null)
  ifndef ROOTCONFIG
    define err_message
Could not find a working ROOT installation.
Please ensure that the root-config executable is on your PATH and try again.
    endef
    $(error "$(err_message)")
  endif
endif

# Figure out the correct shared library file suffix for the current OS
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  SHARED_LIB_SUFFIX=dylib
else ifeq ($(UNAME_S),Linux)
  SHARED_LIB_SUFFIX=so
else
  $(warning Unrecognized operating system encountered.)
  SHARED_LIB_SUFFIX=so
endif

BIN_DIR := bin
INCLUDE_DIR := include
LIB_DIR := lib

ROOT_DICTIONARY := $(LIB_DIR)/dictionaries.o
SHARED_LIB := $(LIB_DIR)/libXSecAnalyzer.$(SHARED_LIB_SUFFIX)

# -MMD -MP make the compiler emit a .d file listing every header each object
# depends on. Without this, editing one of the many header-only implementations
# (CrossSectionExtractor.hh, MCC9SystematicsCalculator.hh, SliceHistogram.hh,
# ...) rebuilds nothing and you silently rerun stale binaries.
CXXFLAGS := $(shell root-config --cflags) -I$(INCLUDE_DIR) -MMD -MP
LDFLAGS := $(shell root-config --libs) \
           -lTMVA -lTMVAGui \
           -L$(LIB_DIR) -lXSecAnalyzer\

ifneq ($(MAKECMDGOALS),debug)
  CXXFLAGS += -O3
else
  CXXFLAGS += -O0 -g
endif

# Source files to use when building the main shared library
SHARED_SOURCES := $(wildcard src/binning/*.cxx)
SHARED_SOURCES += $(wildcard src/selections/*.cxx)
SHARED_SOURCES += $(wildcard src/utils/*.cxx)

# Object files that will be included in the shared library
SHARED_OBJECTS := $(SHARED_SOURCES:.cxx=.o)

# Pin the default goal BEFORE including the dependency files below. The '.d'
# files contain rules of their own, so the first one included would otherwise
# become the default goal and a bare "make" would build a single object file
# instead of the project -- silently leaving every binary stale, which is
# indistinguishable from a successful no-op build.
.DEFAULT_GOAL := all

# Header dependency files emitted by -MMD (one per object, plus one per app
# binary since those compile and link in a single step). Included with a
# leading '-' so a clean tree, where none exist yet, is not an error.
DEPS := $(SHARED_OBJECTS:.o=.d) $(wildcard $(BIN_DIR)/*.d)
-include $(DEPS)

# Causes GNU make to auto-delete the ROOT dictionary object file when the build
# is complete
.INTERMEDIATE: $(ROOT_DICTIONARY)

all: $(SHARED_LIB) bin/ProcessNTuples bin/univmake bin/SlicePlots \
    bin/Unfolder bin/BinScheme bin/StandaloneUnfold bin/xsroot bin/xsnotebook \
    bin/AddFakeWeights bin/AddBeamlineGeometryWeights bin/UnfolderNuMI

debug: all

$(ROOT_DICTIONARY):
	rootcling -f $(LIB_DIR)/dictionaries.cc -c LinkDef.hh
	$(CXX) $(CXXFLAGS) -fPIC -o $@ -c $(LIB_DIR)/dictionaries.cc
	$(RM) $(LIB_DIR)/dictionaries.cc

$(SHARED_OBJECTS): %.o : %.cxx
	$(CXX) $(CXXFLAGS) -fPIC -o $@ -c $<

# -lTMVA is required here, not just in LDFLAGS: CC1mu1piXp uses TMVA::Reader, so
# without it the library carries unresolved TMVA symbols and no DT_NEEDED entry.
# Executables linked with LDFLAGS happened to work, but gSystem->Load() of this
# library (xsroot, xsnotebook, ROOT macros) failed.
$(SHARED_LIB): $(ROOT_DICTIONARY) $(SHARED_OBJECTS)
	$(CXX) $(CXXFLAGS) $(shell root-config --libs) -lTMVA -lTMVAGui \
	  -o $@ -fPIC -shared $^

bin/ProcessNTuples: src/app/ProcessNTuples.C $(SHARED_LIB)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

bin/univmake: src/app/univmake.C $(SHARED_LIB)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

bin/SlicePlots: src/app/Slice_Plots.C $(SHARED_LIB)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

bin/Unfolder: src/app/Unfolder.C $(SHARED_LIB)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

bin/UnfolderNuMI: src/app/UnfolderNuMI.C $(SHARED_LIB)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -O3 -o $@ $<

bin/BinScheme: src/app/binscheme.C $(SHARED_LIB)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

bin/StandaloneUnfold: src/app/standalone_unfold.C $(SHARED_LIB)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

bin/xsroot: $(SHARED_LIB)
	cp src/app/xsroot $(BIN_DIR)

bin/xsnotebook: $(SHARED_LIB)
	cp src/app/xsnotebook $(BIN_DIR)

bin/AddFakeWeights: src/app/NuMI/addFakeWeights.cpp $(SHARED_LIB)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -O3 -o $@ $<

bin/AddBeamlineGeometryWeights: src/app/NuMI/addBeamlineGeometryWeightsToMap.cpp $(SHARED_LIB)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -O3 -o $@ $<

clean:
	$(RM) $(SHARED_LIB) $(BIN_DIR)/*
	$(RM) $(SHARED_OBJECTS)
	$(RM) $(SHARED_OBJECTS:.o=.d)
