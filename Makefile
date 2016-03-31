include Makefile.arch
include Makefile.config

CXX=g++
CXXFLAGS+= -g

CXXFLAGS     += $(ROOTCFLAGS) $(SYSINCLUDES) -O3  -I$(ANITA_UTIL_INSTALL_DIR)/include 
LDFLAGS      += $(ROOTLDFLAGS)  -L$(ANITA_UTIL_INSTALL_DIR)/lib  
LIBS          = $(ROOTLIBS) -g -Wl,-z,defs -lMathMore -lRootFftwWrapper -lAnitaEvent -lAnitaAnalysis
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)
LIBDIR=lib
BUILDDIR=build
INCLUDEDIR=include
BINDIR=bin

.PHONY: clean install all


OBJS := $(addprefix $(BUILDDIR)/, AntennaPositions.o Baseline.o SystemResponse.o UCFilters.o Flags.o Correlator.o dict.o)

#BINARIES := $(addprefix $(BINDIR)/, binary);

INCLUDES := $(addprefix $(INCLUDEDIR)/, $(shell ls $(INCLUDEDIR)))

LINKLIBNAME=UCorrelator
LIBNAME = $(LIBDIR)/lib$(LINKLIBNAME).${DllSuf}

all: $(LIBNAME) $(BINARIES) 

### probably need some magic for Mac OS X here? 
$(LIBNAME): $(OBJS) | $(LIBDIR)
	@echo Building shared library $@
	@$(CXX) $(SOFLAGS) $(LDFLAGS) $(OBJS) $(LIBS) $(GLIBS) -shared -o $@


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


$(BINDIR)/%: %.cc $(INCLUDES) Makefile $(LIBNAME) | $(BINDIR)
	@echo Compiling $<
	@$(CXX)  -I./include -I./ $(CXXFLAGS) -o $@ $(LDFLAGS) -L./$(LIBDIR) -l$(LINKLIBNAME)  $< 

$(BUILDDIR)/dict.cc: $(INCLUDES) LinkDef.h | $(BUILDDIR)
	@echo Running rootcint
	rootcint  -f $@ -c -p -I$(ANITA_UTIL_INSTALL_DIR)/include $(INCLUDES) LinkDef.h

install: 
ifndef ANITA_UTIL_INSTALL_DIR 
	$(error Please define ANITA_UTIL_INSTALL_DIR)
endif 
	install -d $(ANITA_UTIL_INSTALL_DIR)/lib 
	install -d $(ANITA_UTIL_INSTALL_DIR)/cinclude 
	install -c -m 755 $(LIBNAME)(ANITA_UTIL_INSTALL_DIR)/lib  
	install -c -m 644 $(INCLUDES) $(ANITA_UTIL_INSTALL_DIR)/include 

clean: 
	rm -rf build
	rm -rf bin
	rm -rf lib
