###########################################################
#
# Example GNUmakefile to build an Application with Marlin
#
# @author Frank Gaede
# 	
###########################################################

# requires path to MARLIN and Marlin v00-09 or higher
# can be called from top level makefile


#----------------------------------------------------------------------------
#                user variables 
#----------------------------------------------------------------------------


# name of the program/library
PROGNAME = GARLIC


#----  additional user include/libraries
#-- optionally can be specified in the to level GNUmakefile of $MARLIN

USERINCLUDES = -I /home/reinhard/ilcsoft/MarlinUtil/HEAD/include
USERINCLUDES += -I /home/reinhard/ilcsoft/Marlin/HEAD/include
USERINCLUDES += -I ${LCIO}/include
USERINCLUDES += -I ${GEAR}/include -DUSE_GEAR
###USERINCLUDES += -I /home/reinhard/ILCSoft/MarlinReco/Tracking/FullLDCTracking/include

#USERLIBS =
#USERLIBS =  -L /home/reinhard/ILCSoft/MarlinReco/Tracking/FullLDCTracking/lib -lFullLDCTracking
#USERLIBS += -L /home/reinhard/ilcsoft/MarlinUtil/HEAD/lib -lMarlinUtil
#USERLIBS += -Wl,-whole-archive -L /home/reinhard/ILCSoft/MarlinReco/lib -lMarlinReco -Wl,-no-whole-archive

# root stuff
#ROOTLIBS     := $(shell root-config --libs)
#ROOTGLIBS    := $(shell root-config --glibs)

#USERLIBS += $(ROOTLIBS)

#----------------------------------------------------------------------------
#        do not modify anything below in order to use top level makefile
#----------------------------------------------------------------------------

# the base directory of your application
# this is where you  would normaly have the lib, bin and source directories
BASEDIR = ..
LIBDIR = $(BASEDIR)/lib
#BINDIR = $(BASEDIR)/bin

ifdef MARLINDEBUG
 CPPFLAGS = -g -Wall -ansi
else
 CPPFLAGS = -O3 -Wall -ansi
endif


# check for private compiler definition
ifdef MY_CXX
  CXX = $(MY_CXX)
else
  CXX = g++
endif

AR = ar

#--------------------------------------------------------------------------------------
LIBS +=   -Wl,-whole-archive  -L $(LIBDIR)  -l$(PROGNAME) -Wl,-no-whole-archive
#--------------------------------------------------------------------------------------

INCLUDE +=  -I $(BASEDIR)/include
INCLUDE += -I $(ROOTSYS)/include
INCLUDE += $(USERINCLUDES)
CPPFLAGS += $(INCLUDE)
#CPPFLAGS += `. $(MARLIN)/bin/marlin_includes.sh`

#LIBS += $(USERLIBS)
#LIBS += `. $(MARLIN)/bin/marlin_libs.sh`

#SRCS = CopySimHits.cc ECALGarlic.cc ECALPreClustering.cc
objects = $(patsubst %.cc,%.o, $(wildcard *.cc)  )
#objects = $(patsubst %.cc,%.o, $(SRCS)  )

dict = ClusterParametersDict.cc
dictobj = $(patsubst %.cc,%.o, $(dict))

$(LIBDIR)/lib$(PROGNAME).so: $(objects) $(dictobj)
#$(LIBDIR)/lib$(PROGNAME).a: $(objects)
	mkdir -p $(LIBDIR)
	$(CXX) -shared -o $(LIBDIR)/lib$(PROGNAME).so  $?

lib: $(LIBDIR)/lib$(PROGNAME).so

$(dictobj):
	rootcint -f $(dict) -c $(BASEDIR)/include/ClusterParameters.hh-
	$(CXX) -c $(dict) $(CPPFLAGS)

clean:
	rm -f *.o *~ $(LIBDIR)/lib$(PROGNAME).a $(BINDIR)/$(PROGNAME)
	rm -rf $(LIBDIR) $(BINDIR)

doc:
	cd $(BASEDIR)/doc ; doxygen ;

