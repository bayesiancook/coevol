#To compile with gsl, reactivate the following line (remove '#')
USE_GSL=-DGSL

#To statically build the executables, uncomment the following line
#STATIC=-static

ifdef USE_GSL
#GSLLIB=-lgsl -lgslcblas
else
GSLLIB=
endif

#LIBS=$(GSLLIB)
#LIBS=$(GSLLIB) -lboost_regex

CC=g++
#OMP=-fopenmp
#CC=g++
SYSLIB=
INCLUDES=
#CPPFLAGS= $(USE_GSL) $(INCLUDES)
#CPPFLAGS= -g -Wall -O3 $(USE_GSL) $(INCLUDES)
CPPFLAGS= -w -O3 -c -g $(OMP)
#LDFLAGS=$(STATIC)
LDFLAGS= $(OMP)
INSTALL_DIR=
INSTALL_LIB=
SRCS= ConjugateGTRModel.cpp Chrono.cpp BaseType.cpp Chain.cpp DrawTree.cpp GTRSubMatrix.cpp linalg.cpp PhyloProcess.cpp ProbModel.cpp Random.cpp RandomBranchSitePath.cpp RandomTypes.cpp DAGnode.cpp Sample.cpp SequenceAlignment.cpp SiteMapping.cpp BranchSitePath.cpp StateSpace.cpp SubMatrix.cpp TaxonSet.cpp Tree.cpp CodonSubMatrix.cpp MG3OmegaCodonSubMatrix.cpp CodonStateSpace.cpp CodonSequenceAlignment.cpp Move.cpp AminoAcidOmegaSubMatrix.cpp SimilarityMatrix.cpp MGAAProfileMutSelCodonSubMatrix.cpp MGCodonProfileMutSelCodonSubMatrix.cpp MGCodonAAProfileMutSelCodonSubMatrix.cpp BrownianBridge.cpp  PureBrownianProcess.cpp BrownianProcess.cpp BrownianModel.cpp CovMatrix.cpp stabledistribution.cpp

OBJS=$(patsubst %.cpp,%.o,$(SRCS))
ALL_SRCS=$(wildcard *.cpp)
ALL_OBJS=$(patsubst %.cpp,%.o,$(ALL_SRCS))

PROGSDIR=../data
ALL= simplecoevol simplereadcoevol ancov readancov
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

translate$(EXEEXT): $(PROGSDIR)/translate$(EXEEXT)
$(PROGSDIR)/translate$(EXEEXT): Translate.o $(OBJS)
	$(CC) Translate.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

sitemg$(EXEEXT): $(PROGSDIR)/sitemg$(EXEEXT)
$(PROGSDIR)/sitemg$(EXEEXT): SiteMG.o $(OBJS)
	$(CC) SiteMG.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

mixmg$(EXEEXT): $(PROGSDIR)/mixmg$(EXEEXT)
$(PROGSDIR)/mixmg$(EXEEXT): MixMG.o $(OBJS)
	$(CC) MixMG.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

matmixmg$(EXEEXT): $(PROGSDIR)/matmixmg$(EXEEXT)
$(PROGSDIR)/matmixmg$(EXEEXT): MatMixMG.o $(OBJS)
	$(CC) MatMixMG.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

matinfmixmg$(EXEEXT): $(PROGSDIR)/matinfmixmg$(EXEEXT)
$(PROGSDIR)/matinfmixmg$(EXEEXT): MatInfMixMG.o $(OBJS)
	$(CC) MatInfMixMG.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

mutsel$(EXEEXT): $(PROGSDIR)/mutsel$(EXEEXT)
$(PROGSDIR)/mutsel$(EXEEXT): CodonAAProfileMutSelInfMix.o $(OBJS)
	$(CC) CodonAAProfileMutSelInfMix.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

siteom$(EXEEXT): $(PROGSDIR)/siteom$(EXEEXT)
$(PROGSDIR)/siteom$(EXEEXT): SiteOmega.o $(OBJS)
	$(CC) SiteOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

brownian$(EXEEXT): $(PROGSDIR)/brownian$(EXEEXT)
$(PROGSDIR)/brownian$(EXEEXT): Brownian.o $(OBJS)
	$(CC) Brownian.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readbrownian$(EXEEXT): $(PROGSDIR)/readbrownian$(EXEEXT)
$(PROGSDIR)/readbrownian$(EXEEXT): ReadBrownian.o $(OBJS)
	$(CC) ReadBrownian.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readaaprofilemutsel$(EXEEXT): $(PROGSDIR)/readaaprofilemutsel$(EXEEXT)
$(PROGSDIR)/readaaprofilemutsel$(EXEEXT): ReadMutSelMix.o $(OBJS)
	$(CC) ReadMutSelMix.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readaaprofilemutselinfmix$(EXEEXT): $(PROGSDIR)/readaaprofilemutselinfmix$(EXEEXT)
$(PROGSDIR)/readaaprofilemutselinfmix$(EXEEXT): ReadMutSelInfMix.o $(OBJS)
	$(CC) ReadMutSelInfMix.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readcodonaaprofilemutselinfmix$(EXEEXT): $(PROGSDIR)/readcodonaaprofilemutselinfmix$(EXEEXT)
$(PROGSDIR)/readcodonaaprofilemutselinfmix$(EXEEXT): ReadCodonAAProfileMutSelInfMix.o $(OBJS)
	$(CC) ReadCodonAAProfileMutSelInfMix.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readcodonprofilemutselinfmix$(EXEEXT): $(PROGSDIR)/readcodonprofilemutselinfmix$(EXEEXT)
$(PROGSDIR)/readcodonprofilemutselinfmix$(EXEEXT): ReadCodonProfileMutSelInfMix.o $(OBJS)
	$(CC) ReadCodonProfileMutSelInfMix.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readneffmutsel$(EXEEXT): $(PROGSDIR)/readneffmutsel$(EXEEXT)
$(PROGSDIR)/readneffmutsel$(EXEEXT): ReadNeffMutSel.o $(OBJS)
	$(CC) ReadNeffMutSel.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readneffinfmix$(EXEEXT): $(PROGSDIR)/readneffinfmix$(EXEEXT)
$(PROGSDIR)/readneffinfmix$(EXEEXT): ReadNeffInfMix.o $(OBJS)
	$(CC) ReadNeffInfMix.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

neffaaprofileinfmixmutsel$(EXEEXT): $(PROGSDIR)/neffaaprofileinfmixmutsel$(EXEEXT)
$(PROGSDIR)/neffaaprofileinfmixmutsel$(EXEEXT): NeffAAProfileInfMixMutSel.o $(OBJS)
	$(CC) NeffAAProfileInfMixMutSel.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

aaprofilemutselinfmix$(EXEEXT): $(PROGSDIR)/aaprofilemutselinfmix$(EXEEXT)
$(PROGSDIR)/aaprofilemutselinfmix$(EXEEXT): AAProfileMutSelInfMix.o $(OBJS)
	$(CC) AAProfileMutSelInfMix.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

codonprofilemutselinfmix$(EXEEXT): $(PROGSDIR)/codonprofilemutselinfmix$(EXEEXT)
$(PROGSDIR)/codonprofilemutselinfmix$(EXEEXT): CodonProfileMutSelInfMix.o $(OBJS)
	$(CC) CodonProfileMutSelInfMix.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

codonaaprofilemutselinfmix$(EXEEXT): $(PROGSDIR)/codonaaprofilemutselinfmix$(EXEEXT)
$(PROGSDIR)/codonaaprofilemutselinfmix$(EXEEXT): CodonAAProfileMutSelInfMix.o $(OBJS)
	$(CC) CodonAAProfileMutSelInfMix.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

neffaaprofilemutsel$(EXEEXT): $(PROGSDIR)/neffaaprofilemutsel$(EXEEXT)
$(PROGSDIR)/neffaaprofilemutsel$(EXEEXT): NeffAAProfileMutSel.o $(OBJS)
	$(CC) NeffAAProfileMutSel.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

aaprofilemutsel$(EXEEXT): $(PROGSDIR)/aaprofilemutsel$(EXEEXT)
$(PROGSDIR)/aaprofilemutsel$(EXEEXT): AAProfileMutSel.o $(OBJS)
	$(CC) AAProfileMutSel.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

combo$(EXEEXT): $(PROGSDIR)/combo$(EXEEXT)
$(PROGSDIR)/combo$(EXEEXT): CombiCompo.o $(OBJS)
	$(CC) CombiCompo.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readcombo$(EXEEXT): $(PROGSDIR)/readcombo$(EXEEXT)
$(PROGSDIR)/readcombo$(EXEEXT): ReadCombiCompo.o $(OBJS)
	$(CC) ReadCombiCompo.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

redcomp$(EXEEXT): $(PROGSDIR)/redcomp$(EXEEXT)
$(PROGSDIR)/redcomp$(EXEEXT): ReduceCompo.o $(OBJS)
	$(CC) ReduceCompo.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

ic$(EXEEXT): $(PROGSDIR)/ic$(EXEEXT)
$(PROGSDIR)/ic$(EXEEXT): IC.o $(OBJS)
	$(CC) IC.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

compo$(EXEEXT): $(PROGSDIR)/compo$(EXEEXT)
$(PROGSDIR)/compo$(EXEEXT): Compo.o $(OBJS)
	$(CC) Compo.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readcompo$(EXEEXT): $(PROGSDIR)/readcompo$(EXEEXT)
$(PROGSDIR)/readcompo$(EXEEXT): ReadCompo.o $(OBJS)
	$(CC) ReadCompo.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

coevol$(EXEEXT): $(PROGSDIR)/coevol$(EXEEXT)
$(PROGSDIR)/coevol$(EXEEXT): BranchOmegaMultivariate.o $(OBJS)
	$(CC) BranchOmegaMultivariate.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readcoevol$(EXEEXT): $(PROGSDIR)/readcoevol$(EXEEXT)
$(PROGSDIR)/readcoevol$(EXEEXT): ReadBranchOmegaMultivariate.o $(OBJS)
	$(CC) ReadBranchOmegaMultivariate.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

gtr$(EXEEXT): $(PROGSDIR)/gtr$(EXEEXT)
$(PROGSDIR)/gtr$(EXEEXT): GTR.o $(OBJS)
	$(CC) GTR.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

dollo$(EXEEXT): $(PROGSDIR)/dollo$(EXEEXT)
$(PROGSDIR)/dollo$(EXEEXT): Dollo.o $(OBJS)
	$(CC) Dollo.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

simplereadcoevol$(EXEEXT): $(PROGSDIR)/simplereadcoevol$(EXEEXT)
$(PROGSDIR)/simplereadcoevol$(EXEEXT): ReadSimpleBranchOmegaMultivariate.o $(OBJS)
	$(CC) ReadSimpleBranchOmegaMultivariate.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

biclock$(EXEEXT): $(PROGSDIR)/biclock$(EXEEXT)
$(PROGSDIR)/biclock$(EXEEXT): SimpleBiBranchOmegaMultivariate.o $(OBJS)
	$(CC) SimpleBiBranchOmegaMultivariate.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

simplecoevol$(EXEEXT): $(PROGSDIR)/simplecoevol$(EXEEXT)
$(PROGSDIR)/simplecoevol$(EXEEXT): SimpleBranchOmegaMultivariate.o $(OBJS)
	$(CC) SimpleBranchOmegaMultivariate.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

divbd$(EXEEXT): $(PROGSDIR)/divbd$(EXEEXT)
$(PROGSDIR)/divbd$(EXEEXT): DivBD.o $(OBJS)
	$(CC) DivBD.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readdivbd$(EXEEXT): $(PROGSDIR)/readdivbd$(EXEEXT)
$(PROGSDIR)/readdivbd$(EXEEXT): ReadDivBD.o $(OBJS)
	$(CC) ReadDivBD.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

simucont$(EXEEXT): $(PROGSDIR)/simucont$(EXEEXT)
$(PROGSDIR)/simucont$(EXEEXT): SimuCont.o $(OBJS)
	$(CC) SimuCont.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

simu$(EXEEXT): $(PROGSDIR)/simu$(EXEEXT)
$(PROGSDIR)/simu$(EXEEXT): Simu.o $(OBJS)
	$(CC) Simu.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

simucodon$(EXEEXT): $(PROGSDIR)/simucodon$(EXEEXT)
$(PROGSDIR)/simucodon$(EXEEXT): SimuCodon.o $(OBJS)
	$(CC) SimuCodon.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

structstat$(EXEEXT): $(PROGSDIR)/structstat$(EXEEXT)
$(PROGSDIR)/structstat$(EXEEXT): StructStat.o $(OBJS)
	$(CC) StructStat.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

bgc$(EXEEXT): $(PROGSDIR)/bgc$(EXEEXT)
$(PROGSDIR)/bgc$(EXEEXT): BGC.o $(OBJS)
	$(CC) BGC.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

linbgc$(EXEEXT): $(PROGSDIR)/linbgc$(EXEEXT)
$(PROGSDIR)/linbgc$(EXEEXT): LinBGC.o $(OBJS)
	$(CC) LinBGC.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

normalontology$(EXEEXT): $(PROGSDIR)/normalontology$(EXEEXT)
$(PROGSDIR)/normalontology$(EXEEXT): NormalOntology.o $(OBJS)
	$(CC) NormalOntology.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readnormalontology$(EXEEXT): $(PROGSDIR)/readnormalontology$(EXEEXT)
$(PROGSDIR)/readnormalontology$(EXEEXT): ReadNormalOntology.o $(OBJS)
	$(CC) ReadNormalOntology.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

mpiunimixomega$(EXEEXT): $(PROGSDIR)/mpiunimixomega$(EXEEXT)
$(PROGSDIR)/mpiunimixomega$(EXEEXT): UniformMixOmega.o $(OBJS)
	$(CC) UniformMixOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readmpiunimixomega$(EXEEXT): $(PROGSDIR)/readmpiunimixomega$(EXEEXT)
$(PROGSDIR)/readmpiunimixomega$(EXEEXT): ReadUniformMixOmega.o $(OBJS)
	$(CC) ReadUniformMixOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

mpimixomega$(EXEEXT): $(PROGSDIR)/mpimixomega$(EXEEXT)
$(PROGSDIR)/mpimixomega$(EXEEXT): MixOmega.o $(OBJS)
	$(CC) MixOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readmpimixomega$(EXEEXT): $(PROGSDIR)/readmpimixomega$(EXEEXT)
$(PROGSDIR)/readmpimixomega$(EXEEXT): ReadMixOmega.o $(OBJS)
	$(CC) ReadMixOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

mixbgc$(EXEEXT): $(PROGSDIR)/mixbgc$(EXEEXT)
$(PROGSDIR)/mixbgc$(EXEEXT): MixBGC.o $(OBJS)
	$(CC) MixBGC.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readmixbgc$(EXEEXT): $(PROGSDIR)/readmixbgc$(EXEEXT)
$(PROGSDIR)/readmixbgc$(EXEEXT): ReadMixBGC.o $(OBJS)
	$(CC) ReadMixBGC.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readbgc$(EXEEXT): $(PROGSDIR)/readbgc$(EXEEXT)
$(PROGSDIR)/readbgc$(EXEEXT): ReadBGC.o $(OBJS)
	$(CC) ReadBGC.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

life$(EXEEXT): $(PROGSDIR)/life$(EXEEXT)
$(PROGSDIR)/life$(EXEEXT): LifeHistoryRegression.o $(OBJS)
	$(CC) LifeHistoryRegression.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readlife$(EXEEXT): $(PROGSDIR)/readlife$(EXEEXT)
$(PROGSDIR)/readlife$(EXEEXT): ReadLifeHistoryRegression.o $(OBJS)
	$(CC) ReadLifeHistoryRegression.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

bosom$(EXEEXT): $(PROGSDIR)/bosom$(EXEEXT)
$(PROGSDIR)/bosom$(EXEEXT): ConjugateOmega.o $(OBJS)
	$(CC) ConjugateOmega.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

branchmodel$(EXEEXT): $(PROGSDIR)/branchmodel$(EXEEXT)
$(PROGSDIR)/branchmodel$(EXEEXT): Branch.o $(OBJS)
	$(CC) Branch.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readbranchmodel$(EXEEXT): $(PROGSDIR)/readbranchmodel$(EXEEXT)
$(PROGSDIR)/readbranchmodel$(EXEEXT): ReadBranch.o $(OBJS)
	$(CC) ReadBranch.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

fourfold$(EXEEXT): $(PROGSDIR)/fourfold$(EXEEXT)
$(PROGSDIR)/fourfold$(EXEEXT): FourFoldDegenerate.o $(OBJS)
	$(CC) FourFoldDegenerate.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

comp$(EXEEXT): $(PROGSDIR)/comp$(EXEEXT)
$(PROGSDIR)/comp$(EXEEXT): ThirdPositionComposition.o $(OBJS)
	$(CC) ThirdPositionComposition.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

splitdata$(EXEEXT): $(PROGSDIR)/splitdata$(EXEEXT)
$(PROGSDIR)/splitdata$(EXEEXT): SplitData.o $(OBJS)
	$(CC) SplitData.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

partpos$(EXEEXT): $(PROGSDIR)/partpos$(EXEEXT)
$(PROGSDIR)/partpos$(EXEEXT): PartitionCodonPosition.o $(OBJS)
	$(CC) PartitionCodonPosition.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

logcont$(EXEEXT): $(PROGSDIR)/logcont$(EXEEXT)
$(PROGSDIR)/logcont$(EXEEXT): LogCont.o $(OBJS)
	$(CC) LogCont.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

timeline$(EXEEXT): $(PROGSDIR)/timeline$(EXEEXT)
$(PROGSDIR)/timeline$(EXEEXT): TimeLineAncestralCovariance.o $(OBJS)
	$(CC) TimeLineAncestralCovariance.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readtimeline$(EXEEXT): $(PROGSDIR)/readtimeline$(EXEEXT)
$(PROGSDIR)/readtimeline$(EXEEXT): ReadTimeLineAncestralCovariance.o $(OBJS)
	$(CC) ReadTimeLineAncestralCovariance.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

cont$(EXEEXT): $(PROGSDIR)/cont$(EXEEXT)
$(PROGSDIR)/cont$(EXEEXT): Cont.o $(OBJS)
	$(CC) Cont.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readcont$(EXEEXT): $(PROGSDIR)/readcont$(EXEEXT)
$(PROGSDIR)/readcont$(EXEEXT): ReadCont.o $(OBJS)
	$(CC) ReadCont.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

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

tip$(EXEEXT): $(PROGSDIR)/tip$(EXEEXT)
$(PROGSDIR)/tip$(EXEEXT): TipDating.o $(OBJS)
	$(CC) TipDating.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readtip$(EXEEXT): $(PROGSDIR)/readtip$(EXEEXT)
$(PROGSDIR)/readtip$(EXEEXT): ReadTipDating.o $(OBJS)
	$(CC) ReadTipDating.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

logn$(EXEEXT): $(PROGSDIR)/logn$(EXEEXT)
$(PROGSDIR)/logn$(EXEEXT): LogN.o $(OBJS)
	$(CC) LogN.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readlogn$(EXEEXT): $(PROGSDIR)/readlogn$(EXEEXT)
$(PROGSDIR)/readlogn$(EXEEXT): ReadLogN.o $(OBJS)
	$(CC) ReadLogN.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

bdlogn$(EXEEXT): $(PROGSDIR)/bdlogn$(EXEEXT)
$(PROGSDIR)/bdlogn$(EXEEXT): BDLogN.o $(OBJS)
	$(CC) BDLogN.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

drawtree$(EXEEXT): $(PROGSDIR)/drawtree$(EXEEXT)
$(PROGSDIR)/drawtree$(EXEEXT): DrawTreeMain.o $(OBJS)
	$(CC) DrawTreeMain.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

drawcomptree$(EXEEXT): $(PROGSDIR)/drawcomptree$(EXEEXT)
$(PROGSDIR)/drawcomptree$(EXEEXT): DrawCompTree.o $(OBJS)
	$(CC) DrawCompTree.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

simucoevol$(EXEEXT): $(PROGSDIR)/simucoevol$(EXEEXT)
$(PROGSDIR)/simucoevol$(EXEEXT): SimuBranchOmegaMultivariate.o $(OBJS)
	$(CC) SimuBranchOmegaMultivariate.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@


clean:
	-rm -f *.o *.d *.d.*
	-rm -f $(PROGS)

