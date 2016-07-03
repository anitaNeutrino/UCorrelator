include Makefile.arch
include Makefile.config

CXX=g++
CXXFLAGS+= -g

CXXFLAGS     += $(ROOTCFLAGS) $(SYSINCLUDES) -I$(ANITA_UTIL_INSTALL_DIR)/include -I$(EIGEN3_INCLUDE_DIR)/eigen3 
LDFLAGS      += $(ROOTLDFLAGS)  -L$(ANITA_UTIL_INSTALL_DIR)/lib 
LIBS          = $(ROOTLIBS) -g -Wl,-z,defs -lMathMore -lRootFftwWrapper -lAnitaEvent -lAnitaAnalysis -lAnitaCorrelator  
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)
LIBDIR=lib
BUILDDIR=build
INCLUDEDIR=include
BINDIR=bin

.PHONY: clean install all doc


OBJS := $(addprefix $(BUILDDIR)/, AntennaPositions.o Baseline.o SystemResponse.o UCFilters.o Flags.o Correlator.o WaveformCombiner.o PeakFinder.o Analyzer.o AnalysisConfig.o Util.o UCorrelatorDict.o)

BINARIES := $(addprefix $(BINDIR)/, doWais doDecimated doDecimatedNoFilter doLDB);

INCLUDES := $(addprefix $(INCLUDEDIR)/, $(shell ls $(INCLUDEDIR)))

LINKLIBNAME=UCorrelator
LIBNAME = $(LIBDIR)/lib$(LINKLIBNAME).${DllSuf}

all: $(LIBNAME) $(BINARIES) 

### probably need some magic for Mac OS X here? 
$(LIBNAME): $(OBJS) | $(LIBDIR)
	@echo Building shared library $@
	@$(CXX) $(SOFLAGS) $(LDFLAGS) $(OBJS) $(LIBS) $(GLIBS) -shared -o $@
	cp $(BUILDDIR)/*.pcm $(LIBDIR) 


$(OBJS): | $(BUILDDIR)

$(BUILDDIR): 
	mkdir -p $(BUILDDIR)

$(BINDIR): 
	mkdir -p $(BINDIR)

$(LIBDIR): 
	mkdir -p $(LIBDIR)


$(BUILDDIR)/%.o: src/%.cc $(INCLUDES) Makefile | $(BUILDDIR) $(VECTORIZE)
	@echo Compiling  $< 
	@$(CXX)  -I./include $(CXXFLAGS) -o $@ -c $< 

$(BUILDDIR)/%.o: build/%.cc $(INCLUDES) Makefile | $(BUILDDIR) 
	@echo Compiling  $< 
	$(CXX)  -I../include -I./ $(CXXFLAGS) -o $@ -c $< 



$(BINDIR)/%: drivers/%.cc $(INCLUDES) Makefile $(LIBNAME) | $(BINDIR)
	@echo Compiling $<
	@$(CXX)  -I./include -I./ $(CXXFLAGS) -o $@ -L./$(LIBDIR) $(LDFLAGS) -l$(LINKLIBNAME)  $< 

$(BUILDDIR)/UCorrelatorDict.cc: $(INCLUDES) LinkDef.h | $(BUILDDIR)
	@echo Running rootcint
	rootcint  -f $@ -c -p -I$(ANITA_UTIL_INSTALL_DIR)/include $(INCLUDES) LinkDef.h


install: $(LIBNAME)
ifndef ANITA_UTIL_INSTALL_DIR 
	$(error Please define ANITA_UTIL_INSTALL_DIR)
endif 
	install -d $(ANITA_UTIL_INSTALL_DIR)/lib 
	install -d $(ANITA_UTIL_INSTALL_DIR)/include 
	install -c -m 755 $(LIBNAME) $(ANITA_UTIL_INSTALL_DIR)/lib
	install -c -m 644 $(INCLUDES) $(ANITA_UTIL_INSTALL_DIR)/include 
	if [ -e $(BUILDDIR)/UCorrelatorDict_rdict.pcm ];  then install -c -m 755 $(BUILDDIR)/UCorrelatorDict_rdict.pcm $(ANITA_UTIL_INSTALL_DIR)/lib; fi; 



doc: $(INCLUDES) 
	doxygen doc/Doxyfile 
	make -C doc/latex

clean: 
	rm -rf build
	rm -rf bin
	rm -rf lib
