CC=mpic++
SYSLIB=
INCLUDES=
CPPFLAGS= -w -O3 --std=c++11 
#LDFLAGS=$(STATIC)
LDFLAGS= $(OMP)
INSTALL_DIR=
INSTALL_LIB=
SRCS= Chrono.cpp BaseType.cpp Chain.cpp DrawTree.cpp GTRSubMatrix.cpp linalg.cpp PhyloProcess.cpp ProbModel.cpp Random.cpp RandomBranchSitePath.cpp RandomTypes.cpp DAGnode.cpp Sample.cpp SequenceAlignment.cpp SiteMapping.cpp BranchSitePath.cpp StateSpace.cpp SubMatrix.cpp TaxonSet.cpp Tree.cpp CodonSubMatrix.cpp MSCodonSubMatrix.cpp MG3OmegaCodonSubMatrix.cpp CodonStateSpace.cpp CodonSequenceAlignment.cpp Move.cpp AminoAcidOmegaSubMatrix.cpp SimilarityMatrix.cpp CovMatrix.cpp 

OBJS=$(patsubst %.cpp,%.o,$(SRCS))
ALL_SRCS=$(wildcard *.cpp)
ALL_OBJS=$(patsubst %.cpp,%.o,$(ALL_SRCS))

PROGSDIR=../data
ALL= drawtree coevol readcoevol nearlyneutral readnearlyneutral 
PROGS=$(addprefix $(PROGSDIR)/, $(ALL))

# If we are on a windows platform, executables are .exe files
# [S.Bigaret] we assume here that we are compiling w/ MinGW, not Cygwin
EXEEXT=
uname := $(shell uname -s 2>/dev/null)
ifeq (MINGW,$(findstring MINGW,$(uname)))
	EXEEXT=.exe
	PROGS=$(addsuffix $(EXEEXT), $(ALL))
endif

.PHONY: all clean
all: $(PROGS)

# Rules for generate the dependencies automatically

%.d: %.cpp
	@echo "Generating dependencies for $<..."; \
	 set -e; rm -f $@; $(CC) -MM $(CPPFLAGS) $< > $@.$$$$; \
	 sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; rm -f $@.$$$$


# Rules to create .o files from .cpp files
%.o: %.cpp %.d
	$(CC) -c $(CPPFLAGS) $*.cpp

# Include the dependencies unless the request was to clean eveything up
ifneq ($(MAKECMDGOALS),clean)
-include $(ALL_OBJS:.o=.d)
endif

# Targets

coevol$(EXEEXT): $(PROGSDIR)/coevol$(EXEEXT)
$(PROGSDIR)/coevol$(EXEEXT): BranchOmegaMultivariate.o $(OBJS)
	$(CC) BranchOmegaMultivariate.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readcoevol$(EXEEXT): $(PROGSDIR)/readcoevol$(EXEEXT)
$(PROGSDIR)/readcoevol$(EXEEXT): ReadBranchOmegaMultivariate.o $(OBJS)
	$(CC) ReadBranchOmegaMultivariate.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

nearlyneutral$(EXEEXT): $(PROGSDIR)/nearlyneutral$(EXEEXT)
$(PROGSDIR)/nearlyneutral$(EXEEXT): NearlyNeutralModel2.o $(OBJS)
	$(CC) NearlyNeutralModel2.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readnearlyneutral$(EXEEXT): $(PROGSDIR)/readnearlyneutral$(EXEEXT)
$(PROGSDIR)/readnearlyneutral$(EXEEXT): ReadNearlyNeutralModel2.o $(OBJS)
	$(CC) ReadNearlyNeutralModel2.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

$(PROGSDIR)/tipcoevol$(EXEEXT): TipCoevol.o $(OBJS)
	$(CC) TipCoevol.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readtipcoevol$(EXEEXT): $(PROGSDIR)/readtipcoevol$(EXEEXT)
$(PROGSDIR)/readtipcoevol$(EXEEXT): ReadTipCoevol.o $(OBJS)
	$(CC) ReadTipCoevol.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

drawtree$(EXEEXT): $(PROGSDIR)/drawtree$(EXEEXT)
$(PROGSDIR)/drawtree$(EXEEXT): DrawTreeMain.o $(OBJS)
	$(CC) DrawTreeMain.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

clean:
	-rm -f *.o *.d *.d.*
	-rm -f $(PROGS)

