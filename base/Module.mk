# Module.mk for base module
# Copyright (c) 2000 Rene Brun and Fons Rademakers
#
# Author: Fons Rademakers, 29/2/2000

MODDIR       := base
MODDIRS      := $(MODDIR)/src
MODDIRI      := $(MODDIR)/inc

BASEDIR      := $(MODDIR)
BASEDIRS     := $(BASEDIR)/src
BASEDIRI     := $(BASEDIR)/inc

##### libBase (part of libCore) #####
BASEL1       := $(MODDIRI)/LinkDef1.h
BASEL2       := $(MODDIRI)/LinkDef2.h
BASEL3       := $(MODDIRI)/LinkDef3.h
BASEDS1      := $(MODDIRS)/G__Base1.cxx
BASEDS2      := $(MODDIRS)/G__Base2.cxx
BASEDS3      := $(MODDIRS)/G__Base3.cxx
BASEDO1      := $(BASEDS1:.cxx=.o)
BASEDO2      := $(BASEDS2:.cxx=.o)
BASEDO3      := $(BASEDS3:.cxx=.o)
BASEDS       := $(BASEDS1) $(BASEDS2) $(BASEDS3)
BASEDO       := $(BASEDO1) $(BASEDO2) $(BASEDO3)
BASEDH       := $(BASEDS:.cxx=.h)

BASEH1       := $(wildcard $(MODDIRI)/T*.h)
BASEH3       := GuiTypes.h KeySymbols.h
BASEH3       := $(patsubst %,$(MODDIRI)/%,$(BASEH3))
BASEH        := $(filter-out $(MODDIRI)/LinkDef%,$(wildcard $(MODDIRI)/*.h))
BASES        := $(filter-out $(MODDIRS)/G__%,$(wildcard $(MODDIRS)/*.cxx))
BASEO        := $(BASES:.cxx=.o)

BASEDEP      := $(BASEO:.o=.d) $(BASEDO:.o=.d)

# used in the main Makefile
ALLHDRS     += $(patsubst $(MODDIRI)/%.h,include/%.h,$(BASEH))

# include all dependency files
INCLUDEFILES += $(BASEDEP)

##### local rules #####
include/%.h:    $(BASEDIRI)/%.h
		cp $< $@

$(BASEDS1):     $(BASEH1) $(BASEL1) $(ROOTCINTTMP)
		@echo "Generating dictionary $@..."
		$(ROOTCINTTMP) -f $@ -c $(BASEH1) $(BASEL1)
$(BASEDS2):     $(BASEH1) $(BASEL2) $(ROOTCINTTMP)
		@echo "Generating dictionary $@..."
		$(ROOTCINTTMP) -f $@ -c $(BASEH1) $(BASEL2)
$(BASEDS3):     $(BASEH3) $(BASEL3) $(ROOTCINTTMP)
		@echo "Generating dictionary $@..."
		$(ROOTCINTTMP) -f $@ -c $(BASEH3) $(BASEL3)

$(BASEDO1):     $(BASEDS1)
		$(CXX) $(NOOPT) $(CXXFLAGS) -I. -o $@ -c $<
$(BASEDO2):     $(BASEDS2)
		$(CXX) $(NOOPT) $(CXXFLAGS) -I. -o $@ -c $<
$(BASEDO3):     $(BASEDS3)
		$(CXX) $(NOOPT) $(CXXFLAGS) -I. -o $@ -c $<

all-base:       $(BASEO) $(BASEDO)

clean-base:
		@rm -f $(BASEO) $(BASEDO)

clean::         clean-base

distclean-base: clean-base
		@rm -f $(BASEDEP) $(BASEDS) $(BASEDH)

distclean::     distclean-base
