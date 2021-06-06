CC=g++
SYSLIB=
INCLUDES=
CPPFLAGS= -w -O3 -std=c++11
LDFLAGS= 
INSTALL_DIR=
INSTALL_LIB=
SRCS= Chrono.cpp BaseType.cpp Chain.cpp DrawTree.cpp GTRSubMatrix.cpp linalg.cpp PhyloProcess.cpp ProbModel.cpp Random.cpp RandomBranchSitePath.cpp RandomTypes.cpp DAGnode.cpp Sample.cpp SequenceAlignment.cpp SiteMapping.cpp BranchSitePath.cpp StateSpace.cpp SubMatrix.cpp TaxonSet.cpp Tree.cpp CodonSubMatrix.cpp MG3OmegaCodonSubMatrix.cpp CodonStateSpace.cpp CodonSequenceAlignment.cpp Move.cpp AminoAcidOmegaSubMatrix.cpp SimilarityMatrix.cpp BrownianBridge.cpp  PureBrownianProcess.cpp BrownianProcess.cpp BrownianModel.cpp CovMatrix.cpp correlation.cpp correl.cpp

OBJS=$(patsubst %.cpp,%.o,$(SRCS))
ALL_SRCS=$(wildcard *.cpp)
ALL_OBJS=$(patsubst %.cpp,%.o,$(ALL_SRCS))

PROGSDIR=../data
ALL= coevol readcoevol tipcoevol readtipcoevol ancov readancov brownian readbrownian tracecomp drawtree 
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

brownian$(EXEEXT): $(PROGSDIR)/brownian$(EXEEXT)
$(PROGSDIR)/brownian$(EXEEXT): Brownian.o $(OBJS)
	$(CC) Brownian.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readbrownian$(EXEEXT): $(PROGSDIR)/readbrownian$(EXEEXT)
$(PROGSDIR)/readbrownian$(EXEEXT): ReadBrownian.o $(OBJS)
	$(CC) ReadBrownian.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

coevol$(EXEEXT): $(PROGSDIR)/coevol$(EXEEXT)
$(PROGSDIR)/coevol$(EXEEXT): BranchOmegaMultivariate.o $(OBJS)
	$(CC) BranchOmegaMultivariate.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readcoevol$(EXEEXT): $(PROGSDIR)/readcoevol$(EXEEXT)
$(PROGSDIR)/readcoevol$(EXEEXT): ReadBranchOmegaMultivariate.o $(OBJS)
	$(CC) ReadBranchOmegaMultivariate.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

ancov$(EXEEXT): $(PROGSDIR)/ancov$(EXEEXT)
$(PROGSDIR)/ancov$(EXEEXT): AncestralCovariance.o $(OBJS)
	$(CC) AncestralCovariance.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readancov$(EXEEXT): $(PROGSDIR)/readancov$(EXEEXT)
$(PROGSDIR)/readancov$(EXEEXT): ReadAncestralCovariance.o $(OBJS)
	$(CC) ReadAncestralCovariance.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

tipcoevol$(EXEEXT): $(PROGSDIR)/tipcoevol$(EXEEXT)
$(PROGSDIR)/tipcoevol$(EXEEXT): TipCoevol.o $(OBJS)
	$(CC) TipCoevol.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readtipcoevol$(EXEEXT): $(PROGSDIR)/readtipcoevol$(EXEEXT)
$(PROGSDIR)/readtipcoevol$(EXEEXT): ReadTipCoevol.o $(OBJS)
	$(CC) ReadTipCoevol.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

drawtree$(EXEEXT): $(PROGSDIR)/drawtree$(EXEEXT)
$(PROGSDIR)/drawtree$(EXEEXT): DrawTreeMain.o $(OBJS)
	$(CC) DrawTreeMain.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

tracecomp$(EXEEXT): $(PROGSDIR)/tracecomp$(EXEEXT)
$(PROGSDIR)/tracecomp$(EXEEXT): TraceComp.o $(OBJS)
	$(CC) TraceComp.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

clean:
	-rm -f *.o *.d *.d.*
	-rm -f $(PROGS)

